import ropf.Optimization as Optim
import cvxpy as cvx


def test_init_var(init_model):
    cvxvar = Optim.init_var(init_model)
    assert isinstance(cvxvar, dict)
    assert isinstance(cvxvar['U'], cvx.Variable)
    assert isinstance(cvxvar['I'], cvx.Variable)
    assert isinstance(cvxvar['P'], cvx.Variable)
    assert isinstance(cvxvar['Q'], cvx.Variable)
    assert isinstance(cvxvar['Ppv'], cvx.Variable)
    assert isinstance(cvxvar['Pst'], cvx.Variable)
    assert isinstance(cvxvar['Qpv'], cvx.Variable)
    assert isinstance(cvxvar['Qst'], cvx.Variable)
    assert cvxvar['U'].size == (10, 1)
    assert cvxvar['I'].size == (9, 1)
    assert cvxvar['P'].size == (9, 1)
    assert cvxvar['Q'].size == (9, 1)
    assert cvxvar['Ppv'].size == (10, 1)
    assert cvxvar['Pst'].size == (10, 1)
    assert cvxvar['Qpv'].size == (10, 1)
    assert cvxvar['Qst'].size == (10, 1)


def test_init_par(init_model):
    cvxpar = Optim.init_par(init_model)
    assert isinstance(cvxpar, dict)
    assert isinstance(cvxpar['active_load'], cvx.Parameter)
    assert isinstance(cvxpar['reactive_load'], cvx.Parameter)
    assert isinstance(cvxpar['pv_set_points'], cvx.Parameter)
    assert isinstance(cvxpar['storage_set_points'], cvx.Parameter)
    assert isinstance(cvxpar['storage_upper_bounds'], cvx.Parameter)
    assert isinstance(cvxpar['storage_lower_bounds'], cvx.Parameter)
    assert cvxpar['active_load'].size == (10, 1)
    assert cvxpar['reactive_load'].size == (10, 1)
    assert cvxpar['pv_set_points'].size == (10, 1)
    assert cvxpar['storage_set_points'].size == (10, 1)
    assert cvxpar['storage_upper_bounds'].size == (10, 1)
    assert cvxpar['storage_lower_bounds'].size == (10, 1)


def test_power_flow_constraints(init_model):
    cvxpar = Optim.init_par(init_model)
    cvxvar = Optim.init_var(init_model)
    constraints = Optim.power_flow_constraints(init_model, cvxvar, cvxpar)

    for i in constraints:
        isinstance(i, cvx.constraints.eq_constraint.EqConstraint)

    assert len(constraints) == 4*len(init_model.NetworkModel.resistance)


def test_network_limits_constraints(init_model):
    cvxpar = Optim.init_par(init_model)
    cvxvar = Optim.init_var(init_model)
    constraints = Optim.network_limits_constraints(init_model, cvxvar, cvxpar)

    for i in constraints:
        isinstance(i, cvx.constraints.eq_constraint.EqConstraint)

    assert len(constraints) == 3*len(init_model.NetworkModel.resistance)+2


def test_der_injection_constraints(init_model):
    cvxpar = Optim.init_par(init_model)
    cvxvar = Optim.init_var(init_model)
    constraints = Optim.der_injection_constraints(init_model, cvxvar, cvxpar)

    for i in constraints:
        isinstance(i, cvx.constraints.eq_constraint.EqConstraint)

    assert len(constraints) == 6*len(init_model.NetworkModel.resistance)
