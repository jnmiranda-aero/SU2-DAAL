import numpy as np 
import SU2
import SU2.io
import os

# config = SU2.io.Config('config_template.cfg')
# config['MACH_NUMBER'] = 0.69
# config.write('pyWrapperTestConfig.cfg')

# SU2.io.config.dump_config("testConfig.cfg",config)

config = SU2.io.Config()
config = SU2.io.config.dump_config("t.cfg", config)
SU2.io.Config.dump(config,"t.cfg")

# .dump("config.cfg")

# config = SU2.io.Config()