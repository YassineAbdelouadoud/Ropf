# -*- coding: utf-8 -*-
"""

These are the definitions of classes to obtain a complete simulation model

"""

import numpy


class LoadModel:
    """This class contains the description of the loads in the network

    Attributes:
        | active_load (ndarray): An array of active loads, with a column for each bus and a line for each time step.
        | reactive_load (ndarray): An array of reactive loads, with a column for each bus and a line for each time step.
    """

    def __init__(self, input_folder):
        """
        Instance attributes are initialized to None and then method init_loadmodel is called to
        populate the attributes

        Args:
            input_folder (string) : the path to the folder containing the input data
        """

        self.active_load = None
        self.reactive_load = None
        self.init_loadmodel(input_folder)

    def init_loadmodel(self, input_folder):

        """
        This method will populate a LoadModel object attributes by obtaining the relevant data in the input_folder

        Args:
            input_folder (string) : the path to the folder containing the input data

        Returns:
            A LoadModel Object with attributes corresponding to the input data
        """

        input_file = input_folder + "/LoadModel/active_load.csv"
        try:
            active_load = numpy.genfromtxt(input_file, delimiter=',')
            self.active_load = active_load
        except IOError:
            print('The input file active_load.csv could not be located')

        input_file = input_folder + "/LoadModel/reactive_load.csv"
        try:
            reactive_load = numpy.genfromtxt(input_file, delimiter=',')
            self.reactive_load = reactive_load
        except IOError:
            print('The input file active_load.csv could not be located')


class DERModel:
    """This class contains the description of the DERs in the network

    Attributes:
        | pv_inverter_power (ndarray): A vector of installed photovoltaic inverter powers.
        | storage_inverter_power (ndarray): A vector of installed storage inverter powers.
        | pv_set_points (ndarray): An array of photovoltaic set points.
        | storage_set_points (ndarray): A vector of storage set points.
    """

    def __init__(self, input_folder):

        """
        Instance attributes are initialized to None and then method init_loadmodel is called to populate the attributes

        Args:
            input_folder (string) : the path to the folder containing the input data
        """

        self.pv_inverter_power = None
        self.storage_inverter_power = None
        self.pv_set_points = None
        self.storage_set_points = None

        self.init_dermodel(input_folder)

    def init_dermodel(self, input_folder):

        """
        This method will populate a DERModel object attributes by obtaining the relevant data in the input_folder

        Args:
            input_folder (string) : the path to the folder containing the input data

        Returns:
            A DERModel Object with attributes corresponding to the input data
        """

        input_file = input_folder + '/DERModel/pv_inverter.csv'
        try:
            pv_inverter_power = numpy.genfromtxt(input_file, delimiter=',')
            self.pv_inverter_power = pv_inverter_power
        except IOError:
            print('The input file pv_inverter.csv could not be located')

        input_file = input_folder + '/DERModel/storage_inverter.csv'
        try:
            storage_inverter_power = numpy.genfromtxt(input_file, delimiter=',')
            self.storage_inverter_power = storage_inverter_power
        except IOError:
            print('The input file storage_inverter.csv could not be located')

        input_file = input_folder + '/DERModel/pv_set_points.csv'
        try:
            pv_set_points = numpy.genfromtxt(input_file, delimiter=',')
            self.pv_set_points = pv_set_points
        except IOError:
            print('The input file pv_set_points.csv could not be located')

        input_file = input_folder + '/DERModel/storage_set_points.csv'
        try:
            storage_set_points = numpy.genfromtxt(input_file, delimiter=',')
            self.storage_set_points = storage_set_points
        except IOError:
            print('The input file storage_set_points.csv could not be located')


class NetworkModel:
    """This class contains the data relative to the network

        Attributes:
            | from_node (ndarray): A vector of node IDs from which the lines are originating.
            | to_node (ndarray): A vector of node IDs to which the lines are ending.
            | resistance (ndarray): A vector of line resistance.
            | reactance (ndarray): A vector of line reactance.
            | current_limit (ndarray): An array of line current limit.
            | voltage_limit (ndarray): A vector of 2 elements representing the upper and lower bounds on voltage
                magnitude.
            | oltc_constraints (ndarray): A vector of 2 of or more elements. If its length is 2, a continuous model
                will be used for the On-Load Tap Changer, with the two values representing up and down limit of voltage.
                downstream of the OLTC. If it is more than 2, a discrete model will be used, with the values.
                representing the possible voltages
            | downstream of the OLTC

    """

    def __init__(self, input_folder):

        """
        Instance attributes are initialized to None and then method init_networkmodel is called to populate the
        attributes

        Args:
            input_folder (string) : the path to the folder containing the input data
        """

        self.from_node = None
        self.to_node = None
        self.resistance = None
        self.reactance = None
        self.current_limit = None
        self.voltage_limit = None
        self.oltc_constraints = None

        self.init_networkmodel(input_folder)

    def init_networkmodel(self, input_folder):

        """
        This method will populate a NetworkModel object attributes by obtaining the relevant data in the input_folder

        Args:
            input_folder (string) : the path to the folder containing the input data

        Returns:
            A NetworkModel Object with attributes corresponding to the input data
        """

        input_file = input_folder + '/NetworkModel/from_to.csv'
        try:
            from_to = numpy.genfromtxt(input_file, delimiter=',', dtype='int')
            self.from_node = from_to[0]
            self.to_node = from_to[1]
        except IOError:
            print('The input file from_to.csv could not be located')

        input_file = input_folder + '/NetworkModel/reactance.csv'
        try:
            reactance = numpy.genfromtxt(input_file, delimiter=',')
            self.reactance = reactance
        except IOError:
            print('The input file reactance.csv could not be located')

        input_file = input_folder + '/NetworkModel/resistance.csv'
        try:
            resistance = numpy.genfromtxt(input_file, delimiter=',')
            self.resistance = resistance
        except IOError:
            print('The input file resistance.csv could not be located')

        input_file = input_folder + '/NetworkModel/current_limit.csv'
        try:
            current_limit = numpy.genfromtxt(input_file, delimiter=',')
            self.current_limit = current_limit
        except IOError:
            print('The input file current_limit.csv could not be located')

        input_file = input_folder + '/NetworkModel/voltage_limit.csv'
        try:
            voltage_limit = numpy.genfromtxt(input_file, delimiter=',')
            self.voltage_limit = voltage_limit
        except IOError:
            print('The input file voltage_limit.csv could not be located')

        input_file = input_folder + '/NetworkModel/oltc_constraints.csv'
        try:
            oltc_constraints = numpy.genfromtxt(input_file, delimiter=',')
            self.oltc_constraints = oltc_constraints
        except IOError:
            print('The input file oltc_constraints.csv could not be located')


class Model:
    """This class contains the data necessary for a simulation

    Attributes :
        | NetworkModel : A description of the network parameters.
        | LoadModel : A description of the loads in the network.
        | DERModel : A description of the DERs in the network.
    """

    def __init__(self, input_folder):
        self.NetworkModel = NetworkModel(input_folder)
        self.LoadModel = LoadModel(input_folder)
        self.DERModel = DERModel(input_folder)




