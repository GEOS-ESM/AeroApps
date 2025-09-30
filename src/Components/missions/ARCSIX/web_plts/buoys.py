import json
import requests
import csv

access_code = 'X7j3sGXvLFDYuFarvcu0iEPFjXdIy1L1'
 
buoys = {'SIMB3 2024J': '416791b6-2e94-430d-9abf-243ee61928fd',
#         'SIMB3 2024K': 'ad35b741-8057-447f-aeec-9bf77758bbbe',
         'SIMB3 2024L': '60f4bc71-6942-49dc-89ab-d5434424e6dc',
         'SIMB3 2024M': 'b0b375b6-1b79-4ec2-a782-44e13216ec49',
         'SIMB3 2024N': '08e8df7f-cc43-4360-8390-ac2a5afd55df',
         'SIMB3 2024O': '2986a270-4686-4a6c-b608-c28d96bb8e40',
         'SIMB3 2024P': 'a2af67ec-ce23-47e0-98e5-9610248cb605',
         'SIMB3 2024Q': 'b3251af2-037b-4651-a433-ff50e6e30288',
         'SIMB3 2024R': '0ffc82c0-91b4-47fb-a62d-8b16b583196b'
        }

headers = {'Authorization': 'Bearer ' + access_code}
 
url_template = 'https://api.cryosphereinnovation.com/public/deployment/data/{}'

with open('buoys.csv', 'w', newline='') as csvfile:
   spamwriter = csv.writer(csvfile, delimiter=' ',
                           quotechar='|', quoting=csv.QUOTE_MINIMAL)

   for name,id in buoys.items():
 
      url = url_template.format(id,)
      response = requests.get(url, headers=headers).json()
 
      for location in response:
         lon = location['longitude']
         lat = location['latitude']
         time = location['time_stamp']
      spamwriter.writerow([name,time,lon,lat])


