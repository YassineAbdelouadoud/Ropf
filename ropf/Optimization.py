# coding=utf-8
"""
This modules contains the functions used to build the optimization problem and solve it using the iterative procedure
described in `Optimal power flow of a distribution system based on increasingly tight cutting planes added to
a second order cone relaxation <http://www.sciencedirect.com/science/article/pii/S0142061515000095>`__
"""

import cvxpy as cvx
from cvxpy.constraints.second_order import SOC
import numpy


def init_var(model):
    """
    This function starts the set up the optimization problem by defining the variables, which are :
        * U : square of the voltage magnitude
        * I : square of the current
        * P : active power flow
        * Q : reactive power flow
        * Ppv : PV active power injection
        * Qpv : PV reactive power injection
        * Pst : storage active power injection
        * Qst : storage reactive power injection

    Variables relating to nodes have a size equal to the number of nodes in the network, including the root node, even
    if their value at the root node will always be zero (e.g Ppv). Variables relating to line flows have a size equal
    to the number of lines (i.e. number of nodes minus one). Thus, the variables corresponding to the line flow ending
    at node n will be indexed by n-1.

    Args :
        | model : a Model object containing the data necessary for the simulation

    Returns :
        A dictionary containing the variable Objects, with keys being the string corresponding to the variables as
        defined above (U, I, P, Q, Ppv, Pst, Qpv, Qst)
    """

    # Obtaining the number of nodes and lines in the network

    line_cnt = len(model.NetworkModel.resistance)
    node_cnt = len(model.DERModel.pv_inverter_power)

    # Creating the dictionary containing all the optimization variables. Variables are all vectors, with a size
    # depending on the number of lines for variables relating to line flows (I, P and Q) and on the number of nodes
    # for variables corresponding to nodes (U, Ppv, Pst, Qst, Qpv).

    optvar = {'U': cvx.Variable(node_cnt), 'I': cvx.Variable(line_cnt), 'P': cvx.Variable(line_cnt),
              'Q': cvx.Variable(line_cnt), 'Ppv': cvx.Variable(node_cnt), 'Pst': cvx.Variable(node_cnt),
              'Qpv': cvx.Variable(node_cnt), 'Qst': cvx.Variable(node_cnt)}

    return optvar


def init_par(model):
    """
    This function continues the set up the optimization problem by defining the parameters, which are :
        * active_load
        * reactive_load
        * pv_set_points
        * storage_set_points
        * oltc_lim

    Args :
        | model : a Model object containing the data necessary for the simulation

    Returns :
        A dictionary containing the parameter Objects, with keys being the string corresponding to the variables as
        defined above
    """

    # Obtaining the number of nodes in the network

    node_cnt = len(model.DERModel.pv_inverter_power)

    # Creating the dictionary containing all the optimization parameters.

    optpar = {'active_load': cvx.Parameter(node_cnt), 'reactive_load': cvx.Parameter(node_cnt),
              'pv_set_points': cvx.Parameter(node_cnt), 'storage_set_points': cvx.Parameter(node_cnt),
              'oltc_lim': cvx.Parameter(2)}

    return optpar


def set_constraints(model, optvar, optpar):
    """
    This function defines the constraints of the optimization problem

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | optpar : a dictionary containing cvxpy Parameter objects
    Returns :
        A Constraint object
    """

    constraints = power_flow_constraints(model, optvar, optpar) + \
                  network_limits_constraints(model, optvar, optpar) + \
                  der_injection_constraints(model, optvar, optpar)

    return constraints


def power_flow_constraints(model, optvar, optpar):
    """
    This function defines the power flow constraints (i.e. Kirchhoff's laws), corresponding to equations (2a) to (2c)
    along with equation (6) for the relaxed variable change in `Optimal power flow of a distribution system based on
    increasingly tight cutting planes added to a second order cone relaxation
    <http://www.sciencedirect.com/science/article/pii/S0142061515000095/>`__

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | optpar : a dictionary containing cvxpy Parameter objects

    Returns :
        A list of Constraint objects, with a length equal to four times the number of lines
    """

    # Initializing the constraints definition by obtaining the node count and the list of downstream and upstream nodes
    node_cnt = len(model.DERModel.pv_inverter_power)
    upstream_nodes = numpy.concatenate([numpy.zeros(1, dtype='int'), model.NetworkModel.from_node])
    downstream_nodes = []

    for i in range(node_cnt):
        # Obtaining the list of indices for which the current node is in the from vector
        from_idx = numpy.where(model.NetworkModel.from_node == i)
        # Subsetting the previous list of indices into the to_nodes vector to obtain the list of downstream nodes
        downstream_nodes.append(model.NetworkModel.to_node[from_idx])

    constraints = []

    for i in range(1, node_cnt):
        # add the constraints corresponding to the first power flow equation i.e. Kirchhoff's current law on active
        # powers
        constraints.append(optvar['P'][i - 1] - optvar['I'][i - 1] * model.NetworkModel.resistance[i - 1] -
                           sum(optvar['P'][j - 1] for j in downstream_nodes[i]) + optvar['Pst'][i] + optvar['Ppv'][i] ==
                           optpar['active_load'][i])
        # add the constraints corresponding to the second power flow equation i.e. Kirchhoff's current law on reactive
        # powers
        constraints.append(optvar['Q'][i - 1] - optvar['I'][i - 1] * model.NetworkModel.reactance[i - 1] -
                           sum(optvar['Q'][j - 1] for j in downstream_nodes[i]) + optvar['Qst'][i] + optvar['Qpv'][i] ==
                           optpar['reactive_load'][i])
        # add the constraints corresponding to the third power flow equation i.e. Kirchhoff's voltage law
        constraints.append(-optvar['U'][upstream_nodes[i]] + optvar['U'][i] +
                           2 * optvar['P'][i - 1] * model.NetworkModel.resistance[i - 1] +
                           2 * optvar['Q'][i - 1] * model.NetworkModel.reactance[i - 1] -
                           optvar['I'][i - 1] * (model.NetworkModel.resistance[i - 1] ** 2 +
                                                 model.NetworkModel.reactance[i - 1] ** 2) == 0)
        # add the constraint corresponding to the convexified variable change P² + Q² = I*U. The constraint cannot be
        # modelled as it is, because of the I*U product. So we transform the rotated cone to a normal one by rotation
        # and the constraint becomes ||(2*P, 2*Q, U-I)||² - (U+I) <= 0

        constraints.append(cvx.norm(cvx.vstack(2 * optvar['P'][i - 1], 2 * optvar['Q'][i - 1],
                                               optvar['I'][i - 1] - optvar['U'][upstream_nodes[i]])) -
                           optvar['U'][upstream_nodes[i]] - optvar['I'][i - 1] <= 0)

    return constraints


def network_limits_constraints(model, optvar, optpar):
    """
    This function defines the network constraints, corresponding to equations (4a) to (4c)
    in `Optimal power flow of a distribution system based on increasingly tight cutting planes added to a second order
    cone relaxation <http://www.sciencedirect.com/science/article/pii/S0142061515000095/>`__

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | optpar : a dictionary containing cvxpy Parameter objects

    Returns :
        A list of Constraint objects, with a length equal to the number of lines multiplied by three plus two
    """

    # Initializing the constraints definition by obtaining the node count
    node_cnt = len(model.DERModel.pv_inverter_power)

    constraints = []

    for i in range(1, node_cnt):
        # add the constraints corresponding to the limit on line current
        constraints.append(optvar['I'][i - 1] <= model.NetworkModel.current_limit[i - 1] ** 2)
        # add the constraints corresponding to the limits on voltage magnitudes
        constraints.append(optvar['U'][i] <= model.NetworkModel.voltage_limit[1] ** 2)
        constraints.append(optvar['U'][i] >= model.NetworkModel.voltage_limit[0] ** 2)

    constraints.append(optvar['U'][0] <= optpar['oltc_lim'][1] ** 2)
    constraints.append(optvar['U'][0] >= optpar['oltc_lim'][0] ** 2)

    return constraints


def der_injection_constraints(model, optvar, optpar):
    """
    This function defines the constraints on DER injection, corresponding to equations (5a) to (5d)
    in `Optimal power flow of a distribution system based on increasingly tight cutting planes added to a second order
    cone relaxation <http://www.sciencedirect.com/science/article/pii/S0142061515000095/>`__

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | optpar : a dictionary containing cvxpy Parameter objects

    Returns :
        A list of Constraint objects, with a length equal to the number of lines multiplied by 4
    """

    # Initializing the constraints definition by obtaining the node count
    node_cnt = len(model.DERModel.pv_inverter_power)

    constraints = []

    for i in range(1, node_cnt):
        # add the constraints corresponding to the limit on apparent power of PV inverters
        constraints.append(cvx.square(optvar['Ppv'][i]) + cvx.square(optvar['Qpv'][i]) <=
                           model.DERModel.pv_inverter_power[i] ** 2)
        # add the constraints corresponding to the limit on apparent power of storage inverters
        constraints.append(cvx.square(optvar['Pst'][i]) + cvx.square(optvar['Qst'][i]) <=
                           model.DERModel.storage_inverter_power[i] ** 2)
        # add the constraints corresponding to upper bound on PV active injection
        constraints.append(optvar['Ppv'][i] <= optpar['pv_set_points'][i])
        # add the constraints corresponding to lower bound on PV active injection
        constraints.append(optvar['Ppv'][i] >= 0)

    return constraints


def set_objective(model, optvar, optpar, objtype):
    """
    This function defines the objective function, which is a minimization of sum of the losses in the network (r*I) and
    a term representing the distance to the set-points, either with a L-1 or L-2 norm

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | optpar : a dictionary containing cvxpy Parameter objects
        | objtype : a string representing the type of norm used for measuring the distance to the set-points

    Returns :
        A cvxpy Objective object
    """

    if objtype == 'L1':
        obj = cvx.Minimize(cvx.sum_entries(cvx.mul_elemwise(model.NetworkModel.resistance, optvar['I'])) +
                           cvx.norm(optvar['Ppv'] - optpar['pv_set_points'], 1) +
                           cvx.norm(optvar['Pst'] - optpar['storage_set_points'], 1))
    elif objtype == 'L2':
        obj = cvx.Minimize(cvx.sum_entries(cvx.mul_elemwise(model.NetworkModel.resistance, optvar['I'])) +
                           cvx.norm(optvar['Ppv'] - optpar['pv_set_points'], 2) +
                           cvx.norm(optvar['Pst'] - optpar['storage_set_points'], 2))
    elif objtype == 'Loss':
        obj = cvx.Minimize(cvx.sum_entries(cvx.mul_elemwise(model.NetworkModel.resistance, optvar['I'])))

    else:
        raise IOError("Invalid objtype parameter passed to opf_solve, choose between L1 and L2")

    return obj


def parameter_update(optpar, model, step):
    """
    This function updates the value of the parameters for the targeted time step
    :param optpar:
    :param model:
    :param step:
    :return:
    """

    optpar['active_load'].value = model.LoadModel.active_load[step]
    optpar['reactive_load'].value = model.LoadModel.reactive_load[step]
    optpar['pv_set_points'].value = model.DERModel.pv_set_points[step]
    optpar['storage_set_points'].value = model.DERModel.storage_set_points[step]
    optpar['oltc_lim'].value = model.NetworkModel.oltc_constraints

    return optpar


def single_opf_solve(optim):
    """
    This function solves a single step optimal power flow

    Args :
        | optim : a cvxpy Optimization object containing the problem definition

    Returns :
        A list containing the results of the single step optimal power flow
    """


def opf_solve(model, objtype, step=None):
    """
    This function creates an optimization problem cvxpy object by initializing the variable and parameter definitions
    and then calling set_constraints and set_obj

    """

    # initializing the Variable and Parameter cvxpy objects
    optvar = init_var(model)
    optpar = init_par(model)

    # defining the constraints and objective
    constraints = set_constraints(model, optvar, optpar)
    objective = set_objective(model, optvar, optpar, objtype)

    opf_prob = cvx.Problem(objective, constraints)

    optpar = parameter_update(optpar, model, step)

    opf_prob.solve(verbose=True)

    return [opf_prob]