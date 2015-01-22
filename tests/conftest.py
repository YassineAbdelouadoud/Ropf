import os
import pytest
import ropf.Simulation as Sim


@pytest.fixture(scope="session")
def init_model():
    input_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir, "data"))
    model = Sim.Model(input_folder)
    return model