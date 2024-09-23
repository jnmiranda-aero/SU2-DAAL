import os
import numpy as np
import SU2
import SU2.io 

# Generate a blank config file based on the template
config = SU2.io.Config('/home/juan/DAAL/SU2/SU2/config_template.cfg')

SU2.io.config.dump_config("testConfig.cfg",config)

