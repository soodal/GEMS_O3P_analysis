#!/usr/bin/env python
# coding: utf-8

from pywoudc import WoudcClient
import matplotlib.pyplot as plt

import datetime

# get_ipython().run_line_magic('matplotlib', 'inline')


# connect to WOUDC data service
client = WoudcClient()

# get Edmonton sonde data from 1980 to present
# gaw_id = 'EDT'
gaw_id = 'TKB'
begin = datetime.date(1980, 1, 1)
end = datetime.date.today()

bbox = [80, 150, -10, 50]
print(client.get
data = client.get_data('ozonesonde',
                       filters={'gaw_id': gaw_id},
                       temporal=[begin, end])

items = data.keys()

len(data['features'])


type(data['features'])




data['features'][0]


i = -1
csv_content = data['features'][i]['properties']['data_block']
csv_fn = './woudc_ozonesonde_properties_data_block.csv'
with open(csv_fn, 'w') as f:
    f.write(csv_content)
    f.close()
    
df = pd.read_csv(cvs_fn)



# setup graph axes
x_axis = [datetime.datetime.strptime(x['properties']['instance_datetime'], '%Y/%m/%d %H:%M:%S+00')
          for x in data['features'] if x['properties']['flight_summary_totalo3']]
y_axis = [float(x['properties']['flight_summary_totalo3'])
          for x in data['features'] if x['properties']['flight_summary_totalo3']]

# render simple plot
plt.plot(x_axis, y_axis)
plt.savefig('poh.png')
plt.close()



# In[ ]:
