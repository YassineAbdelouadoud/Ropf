# coding=utf-8
"""
This modules contains the functions used to build the optimization problem and solve it using the iterative procedure
described in `Optimal power flow of a distribution system based on increasingly tight cutting planes added to
a second order cone relaxation <http://www.sciencedirect.com/science/article/pii/S0142061515000095>`__
"""

import cvxpy as cvx
import numpy
import powerflow
import ecos


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
        * cut_rhs

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
              'oltc_lim': cvx.Parameter(2), 'cut_rhs': cvx.Parameter(1)}

    return optpar


def set_constraints(model, optvar, optpar, vtol):
    """
    This function defines the constraints of the optimization problem

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | optpar : a dictionary containing cvxpy Parameter objects
        | vtol : a decrease of the upper-bound of voltage magnitude to allow the attainment of a feasible solution
    Returns :
        A Constraint object
    """

    constraints = power_flow_constraints(model, optvar, optpar) + \
                  network_limits_constraints(model, optvar, optpar, vtol) + \
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


def network_limits_constraints(model, optvar, optpar, vtol):
    """
    This function defines the network constraints, corresponding to equations (4a) to (4c)
    in `Optimal power flow of a distribution system based on increasingly tight cutting planes added to a second order
    cone relaxation <http://www.sciencedirect.com/science/article/pii/S0142061515000095/>`__

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | optpar : a dictionary containing cvxpy Parameter objects
        | vtol : a decrease of the upper-bound of voltage magnitude to allow the attainment of a feasible solution

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
        constraints.append(optvar['U'][i] <= (model.NetworkModel.voltage_limit[1] - vtol) ** 2)
        constraints.append(optvar['U'][i] >= (model.NetworkModel.voltage_limit[0] - vtol) ** 2)

    # add the constraints on OLTC voltage magnitude
    constraints.append(optvar['U'][0] <= optpar['oltc_lim'][1] ** 2)
    constraints.append(optvar['U'][0] >= optpar['oltc_lim'][0] ** 2)
    # add the constraint corresponding to the cutting plane sum r*I < sum (P²+Q²)/V²
    constraints.append(cvx.sum_entries(cvx.mul_elemwise(model.NetworkModel.resistance, optvar['I'])) <=
                       optpar['cut_rhs'])

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


def set_objective(model, optvar, optpar):
    """
    This function defines the objective function, which is a minimization of sum of the losses in the network (r*I) and
    a term representing the distance to the set-points, either with a L-1 or L-2 norm

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | optpar : a dictionary containing cvxpy Parameter objects
        | objtype : a string representing the type of objective (can take the values 'L1', 'L2' and 'Loss')

    Returns :
        A cvxpy Objective object
    """

    objtype = 'L1'

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
        raise IOError("Invalid objtype parameter passed to opf_solve, choose between L1, L2 and Loss")

    return obj


def parameter_update(optpar, model, step):
    """
    This function updates the value of the parameters for the targeted time step

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | step : the step chosen for the calculation

    Returns :
        An optpar dictionary containing updated parameter values
    """

    optpar['active_load'].value = model.LoadModel.active_load[step]
    optpar['reactive_load'].value = model.LoadModel.reactive_load[step]
    optpar['pv_set_points'].value = model.DERModel.pv_set_points[step]
    optpar['storage_set_points'].value = model.DERModel.storage_set_points[step]
    optpar['oltc_lim'].value = model.NetworkModel.oltc_constraints
    optpar['cut_rhs'].value = sum(model.NetworkModel.resistance * (model.NetworkModel.current_limit ** 2))

    return optpar


def ecos_data_update(ecos_data, model, step):
    """
    This function updates the value of data in ecos data structure

    Args :
        | ecos_data : optimization data in ecos format
        | model : a Model object containing the data necessary for the simulation
        | step : the step chosen for the calculation

    Returns :
        An optpar dictionary containing updated parameter values
    """

    node_cnt = len(model.DERModel.pv_inverter_power)

    # updating the active power load
    ecos_data[5][0::3] = model.LoadModel.active_load[step,1:]
    ecos_data[5][1::3] = model.LoadModel.reactive_load[step,1:]
    ecos_data[5][2::3] = numpy.zeros(node_cnt - 1)

    # updating the set-points for pv and storage
    ecos_data[2][0:node_cnt] = -model.DERModel.pv_set_points[step]*3
    ecos_data[2][node_cnt:(2 * node_cnt)] = model.DERModel.pv_set_points[step]
    ecos_data[2][2 * node_cnt:(3 * node_cnt)] = -model.DERModel.storage_set_points[step]
    ecos_data[2][3 * node_cnt:(4 * node_cnt)] = model.DERModel.storage_set_points[step]
    ecos_data[2][(8 * node_cnt + 1):(12 * node_cnt):4] = model.DERModel.pv_set_points[step]

    return ecos_data


def get_results(model, optvar, cpp_model, step, optpar):
    """
    This function runs a power flow against the results of the OPF and then organizes the results into a dictionary

    Args :
        | model : a Model object containing the data necessary for the simulation
        | optvar : a dictionary containing cvxpy Variable objects
        | cpp_model : a C++ object containing a network model suitable for powerflow calculations
        | step : time step at which the opf is solved
        | optpar : a dictionary containing cvxpy Parameter objects


    Returns :
        A dictionary containing the results of an optimization problem
    """

    # preparing the data for the power flow simulation
    node_cnt = len(model.DERModel.pv_inverter_power)
    u_opf = numpy.squeeze(numpy.asarray(optvar['U'].value))
    p_opf = numpy.squeeze(numpy.asarray(optvar['P'].value))
    q_opf = numpy.squeeze(numpy.asarray(optvar['Q'].value))
    i_opf = numpy.squeeze(numpy.asarray(optvar['I'].value))
    ppv = numpy.squeeze(numpy.asarray(optvar['Ppv'].value))
    qpv = numpy.squeeze(numpy.asarray(optvar['Qpv'].value))
    pst = numpy.squeeze(numpy.asarray(optvar['Pst'].value))
    qst = numpy.squeeze(numpy.asarray(optvar['Qst'].value))

    u = numpy.ones(node_cnt) * u_opf[0]
    p = model.LoadModel.active_load[step] - ppv - pst
    q = model.LoadModel.reactive_load[step] - qpv - qst

    # launching the BFS powerflow algorithm
    respf = cpp_model.BackwardForward_Sweep(u, p, q, 1e-8)

    # intermediate calculations for some results
    i_opf_real = (p_opf ** 2 + q_opf ** 2) / (u_opf[model.NetworkModel.from_node])
    loss_opf = sum(model.NetworkModel.resistance * i_opf)
    loss_opf_real = sum(model.NetworkModel.resistance * i_opf_real)
    gap = loss_opf - loss_opf_real
    dpv = sum(abs(numpy.squeeze(numpy.asarray(optpar['pv_set_points'].value)) - ppv))
    dst = sum(abs(numpy.squeeze(numpy.asarray(optpar['storage_set_points'].value)) - pst))

    result = {'V_opf': numpy.sqrt(u_opf), 'I_opf': i_opf, 'P_opf': p_opf, 'Q_opf': q_opf,
              'Ppv': ppv, 'Qpv': qpv, 'Pst': pst, 'Qst': qst, 'dPV': dpv, 'dST': dst,
              'V_pf': numpy.sqrt(respf[2]), 'P_pf': respf[0], 'Q_pf': respf[1],
              'Loss_opf': loss_opf, 'Loss_opf_real': loss_opf_real, 'Gap': gap}

    return result


def single_opf_solve(model, opf_prob, ecos_data, optvar, optpar,
                     cpp_model, solver_param, max_iter, t, high_level_verbose):
    """
    This function solves a single step optimal power flow

    Args :
        | model : a Model object containing the data necessary for the simulation
        | opf_prob :
        | ecos_data :
        | optvar : a dictionary containing cvxpy Variable objects
        | optpar : a dictionary containing cvxpy Parameter objects
        | cpp_model : a C++ object containing a network model suitable for powerflow calculations
        | solver_param : a dictionary of solver parameters (defaults to
        |               {'solver': 'ECOS', 'verbose': False, 'tol': 1e-7})
        | max_iter : the maximum number of iterations
        | t : the time step for which the optimization is made
        | high_level_verbose : level of screen output (can be 0, 1 or 2)

    Returns :
        A dictionary containing the results of the single step optimal power flow
    """

    #  we start the iterative optimization
    iter_cnt = 0
    node_cnt = len(model.DERModel.pv_inverter_power)

    # we print the column names of the detailed console output
    if high_level_verbose == 2:
        print '%-5s%-8s%-8s%-8s%-8s%-10s' % ('Iter','dPV','dST','Losses','Gap','Overvoltage')

    # We emulate a do{} while with while True: if test: break
    while True:
        iter_cnt += 1
        solver_output = ecos.solve(ecos_data[0], ecos_data[1], ecos_data[2], ecos_data[3], ecos_data[4], ecos_data[5],
                                   abstol=solver_param['tol'], feastol=solver_param['tol'],
                                   verbose=solver_param['verbose'])
        opf_prob.unpack_results(cvx.ECOS, solver_output)

        # the loop is broken is the status of the solution is either infeasible or unbounded
        if (opf_prob.status == cvx.INFEASIBLE) or (opf_prob.status == cvx.INFEASIBLE_INACCURATE):
            sol_status = 'Infeasible'
            result = dict()
            break

        if (opf_prob.status == cvx.UNBOUNDED) or (opf_prob.status == cvx.UNBOUNDED_INACCURATE):
            sol_status = 'Unbounded'
            result = dict()
            break

        # the results are extracted and formatted
        result = get_results(model, optvar, cpp_model, t, optpar)

        if high_level_verbose == 2:
            print "%-5i%-8.3f%-8.3f%-8.3f%-8.4f%-10.4f" % (iter_cnt, result['dPV'], result['dST'],
                                                           result['Loss_opf'], result['Gap'],
                                                           (max(result['V_pf']) - model.NetworkModel.voltage_limit[1]))
        # the right-hand side of the cut is updated
        ecos_data[2][node_cnt * 8 - 2] = result['Loss_opf_real']

        # the loop is broken if either a feasible solution is attained or the maximum number of iteration is surpassed
        if max(result['V_pf']) < model.NetworkModel.voltage_limit[1]:
            sol_status = 'Optimal'
            break

        if iter_cnt == max_iter:
            sol_status = 'Max iteration attained'
            break

    result['sol_status'] = sol_status
    result['iter_cnt'] = iter_cnt

    return result


def opf_solve(model, step=None, solver_param=None, max_iter=50, vtol=1e-5, high_level_verbose=2):
    """
    This function solves an optimal power flow
    Args :
        | model : a Model object containing the data necessary for the simulation
        | step : a scalar or vector of steps for which the OPF is to be solved
        | max_iter : the limit on SOCP iterations
        | vtol : a decrease of the upper-bound of voltage magnitude to allow the attainment of a feasible solution
        | high_level_verbose : level of screen output (can be 0, 1 or 2)

    Returns :
        A dictionary containing the results of the optimal power flow
    """
    if not solver_param:
        solver_param = {'verbose': False, 'tol': 1e-7}

    # initializing the Variable and Parameter cvxpy objects
    optvar = init_var(model)
    optpar = init_par(model)

    # defining the constraints and objective
    constraints = set_constraints(model, optvar, optpar, vtol)
    objective = set_objective(model, optvar, optpar)

    # setting up the optimization problems
    opf_prob = cvx.Problem(objective, constraints)

    # setting up the cpp object to calculate power flows
    cpp_model = powerflow.PyNetwork(model.NetworkModel.resistance, model.NetworkModel.reactance,
                                    model.NetworkModel.from_node, model.NetworkModel.to_node)

    # setting the value of the parameters for the initial problem

    if isinstance(step, int):
        optpar = parameter_update(optpar, model, step)
    else:
        optpar = parameter_update(optpar, model, step[0])

    # getting the problem data in ecos format
    ecos_data = opf_prob.get_problem_data(cvx.ECOS)

    if isinstance(step, int):
        result = single_opf_solve(model, opf_prob, ecos_data, optvar, optpar,
                                  cpp_model, solver_param, max_iter, step, high_level_verbose)
    else:
        result = dict()
        for t in step:
            ecos_data = ecos_data_update(ecos_data, model, t)
            result[t] = single_opf_solve(model, opf_prob, ecos_data, optvar, optpar,
                                         cpp_model, solver_param, max_iter, t, high_level_verbose)

    return result

