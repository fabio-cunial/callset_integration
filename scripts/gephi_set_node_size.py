# Assuming every node has an attribute 'Size' in the DOT file.
# In Gephi: Window > Console, then:
# execfile('/Users/fcunial/Downloads/gephi_get_node_size.py')
#
from java.awt import Color
import math

for v in g.nodes:
	v.size = int(math.sqrt(float(v.Size)))
