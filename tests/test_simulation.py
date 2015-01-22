import ropf.Simulation as Sim


def test_model_constructor(init_model):
    assert isinstance(init_model, Sim.Model)


def test_dermodel(init_model):
    assert init_model.DERModel.pv_inverter_power.shape == (10, )
    assert init_model.DERModel.storage_inverter_power.shape == (10, )
    assert init_model.DERModel.pv_set_points.shape == (24, 10)
    assert init_model.DERModel.storage_set_points.shape == (24, 10)


def test_loadmodel(init_model):
    assert init_model.LoadModel.active_load.shape == (24, 10)
    assert init_model.LoadModel.reactive_load.shape == (24, 10)


def test_networkmodel(init_model):
    assert init_model.NetworkModel.resistance.shape == (9,)
    assert init_model.NetworkModel.reactance.shape == (9,)
    assert init_model.NetworkModel.from_node.shape == (9,)
    assert init_model.NetworkModel.to_node.shape == (9,)
    assert init_model.NetworkModel.voltage_limit.shape == (2,)
    assert init_model.NetworkModel.oltc_constraints.shape == (2,)