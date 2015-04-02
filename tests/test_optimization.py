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
    assert isinstance(cvxpar['oltc_lim'], cvx.Parameter)
    assert cvxpar['active_load'].size == (10, 1)
    assert cvxpar['reactive_load'].size == (10, 1)
    assert cvxpar['pv_set_points'].size == (10, 1)
    assert cvxpar['storage_set_points'].size == (10, 1)
    assert cvxpar['oltc_lim'].size == (2, 1)


def test_power_flow_constraints(init_model):
    cvxpar = Optim.init_par(init_model)
    cvxvar = Optim.init_var(init_model)
    constraints = Optim.power_flow_constraints(init_model, cvxvar, cvxpar)

    for i in [x for x in range(36) if (x + 1) % 4 != 0]:
        assert isinstance(constraints[i], cvx.constraints.eq_constraint.EqConstraint)

    for i in [x for x in range(36) if (x + 1) % 4 == 0]:
        assert isinstance(constraints[i], cvx.constraints.eq_constraint.LeqConstraint)

    assert len(constraints) == 4 * len(init_model.NetworkModel.resistance)


def test_network_limits_constraints(init_model):
    cvxpar = Optim.init_par(init_model)
    cvxvar = Optim.init_var(init_model)
    constraints = Optim.network_limits_constraints(init_model, cvxvar, cvxpar)

    for i in constraints:
        assert isinstance(i, cvx.constraints.eq_constraint.LeqConstraint)

    assert len(constraints) == 3 * len(init_model.NetworkModel.resistance) + 3


def test_der_injection_constraints(init_model):
    cvxpar = Optim.init_par(init_model)
    cvxvar = Optim.init_var(init_model)
    constraints = Optim.der_injection_constraints(init_model, cvxvar, cvxpar)

    for i in constraints:
        assert isinstance(i, cvx.constraints.eq_constraint.LeqConstraint)

    assert len(constraints) == 4 * len(init_model.NetworkModel.resistance)


def test_parameter_update(init_model):
    cvxpar = Optim.init_par(init_model)
    step = 0
    cvxpar = Optim.parameter_update(cvxpar, init_model, step)

    assert isinstance(cvxpar, dict)
    assert isinstance(cvxpar['active_load'], cvx.Parameter)
    assert isinstance(cvxpar['reactive_load'], cvx.Parameter)
    assert isinstance(cvxpar['pv_set_points'], cvx.Parameter)
    assert isinstance(cvxpar['storage_set_points'], cvx.Parameter)
    assert isinstance(cvxpar['oltc_lim'], cvx.Parameter)
    assert cvxpar['active_load'].size == (10, 1)
    assert cvxpar['reactive_load'].size == (10, 1)
    assert cvxpar['pv_set_points'].size == (10, 1)
    assert cvxpar['storage_set_points'].size == (10, 1)
    assert cvxpar['oltc_lim'].size == (2, 1)


def test_set_objective(init_model):
    cvxvar = Optim.init_var(init_model)
    cvxpar = Optim.init_par(init_model)
    objtype = 'L1'

    obj = Optim.set_objective(init_model, cvxvar, cvxpar, objtype)

    assert isinstance(obj, cvx.problems.objective.Minimize)


def test_opf_solve(init_model):

    solver_param = {'solver': 'ECOS', 'verbose': True, 'tol': 1e-5}
    m = Optim.opf_solve(init_model, objtype='L1', step=0, solver_param=solver_param)