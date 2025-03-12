import numpy as np
from utility import *
class simulation():
    voltage_drive = False
    scan_frequency = False
    mode = None
    lambda_incident = 0
    def __init__(self,lambda_incident):
        
        self.lambda_incident = lambda_incident
        print(self.lambda_incident)
        
    def set_dt(self, *args, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")

    def create_time_array(self, *args, **kwargs):
        raise NotImplementedError("Subclasses must implement this method")