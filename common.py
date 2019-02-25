from __future__ import print_function
from __future__ import division
import sys,os
print('###############################################################################################')
print( "Python:",os.path.realpath(sys.executable), sys.version.split('\n')[0])
PY3 = sys.version_info > (3,)
if(PY3):
  print("\n   !!!   Warning: Python3 detected     !!!")
  
  print("   !!! This program is written python2 !!!\n")



def save_safe_yaml(data,filename):
    yamltext=yaml.safe_dump(data , default_flow_style=False)
    with open(filename, 'w') as outfile:
        yaml.safe_dump(data, outfile, default_flow_style=False)
    return yamltext
def save_yaml(data,filename):
    with open(filename, 'w') as outfile:
        yaml.dump(data, outfile, default_flow_style=False)
