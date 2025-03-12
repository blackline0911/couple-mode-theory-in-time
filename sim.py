import numpy as np
from utility import *
class simulation():
    voltage_drive = False
    scan_frequency = False
    mode = None
    def __init__(self,mode):
        self.mode = mode
        modes = {"voltage_drive", "scan_frequency"}
        assert mode in modes, f"Please choose a simulation mode from {modes}"

    def set_dt(self, *args, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")

    def create_time_array(self, *args, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")