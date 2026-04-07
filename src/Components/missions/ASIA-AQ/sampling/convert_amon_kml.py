#!/usr/bin/env python3
"""
Convert AMoN KML file to CSV format.
"""

import re
import csv
import sys
from xml.etree import ElementTree as ET
from datetime import datetime
from html.parser import HTMLParser
import argparse

class TableParser(HTMLParser):
    """Parse HTML table from KML description field."""
    def __init__(self):
        super().__init__()
        self.data = {}
        self.current_key = None
        self.in_td = False
        self.td_count = 0
        self.current_row_key = None

    def handle_starttag(self, tag, attrs):
        if tag == 'tr':
            self.td_count = 0
            self.current_row_key = None
        if tag == 'td':
            self.in_td = True
            self.td_count += 1

    def handle_endtag(self, tag):
        if tag == 'td':
            self.in_td = False

    def handle_data(self, data):
        data = data.strip()
        if not data or not self.in_td:
            return
        if self.td_count == 1:
            self.current_row_key = data.rstrip(':')
        elif self.td_count == 2 and self.current_row_key:
            self.data[self.current_row_key] = data


def parse_date(date_str):
    """Parse date string like 'Oct 15 2019 12:00AM' to 'YYYY-MM-DD'."""
    if not date_str or not date_str.strip():
        return ''
    date_str = date_str.strip()
    # remove extra spaces
    date_str = re.sub(r'\s+', ' ', date_str)
    try:
        dt = datetime.strptime(date_str, '%b %d %Y %I:%M%p')
        return dt.strftime('%Y-%m-%d')
    except ValueError:
        try:
            dt = datetime.strptime(date_str, '%b  %d %Y %I:%M%p')
            return dt.strftime('%Y-%m-%d')
        except ValueError:
            return date_str


def get_state_info(site_id):
    """
    Derive state abbreviation from site ID prefix.
    Returns (state_abbrev, state_name).
    Note: KML does not contain county or stateName — 
    these are not available in this file.
    """
    # state prefix from site ID (first 2 chars)
    state_map = {
        'AB': ('AB', 'Alberta'),
        'AK': ('AK', 'Alaska'),
        'AL': ('AL', 'Alabama'),
        'AR': ('AR', 'Arkansas'),
        'AZ': ('AZ', 'Arizona'),
        'CA': ('CA', 'California'),
        'CO': ('CO', 'Colorado'),
        'CT': ('CT', 'Connecticut'),
        'FL': ('FL', 'Florida'),
        'GA': ('GA', 'Georgia'),
        'ID': ('ID', 'Idaho'),
        'IL': ('IL', 'Illinois'),
        'IN': ('IN', 'Indiana'),
        'KS': ('KS', 'Kansas'),
        'KY': ('KY', 'Kentucky'),
        'MD': ('MD', 'Maryland'),
        'ME': ('ME', 'Maine'),
        'MI': ('MI', 'Michigan'),
        'MN': ('MN', 'Minnesota'),
        'MS': ('MS', 'Mississippi'),
        'NC': ('NC', 'North Carolina'),
        'NE': ('NE', 'Nebraska'),
        'NH': ('NH', 'New Hampshire'),
        'NJ': ('NJ', 'New Jersey'),
        'NM': ('NM', 'New Mexico'),
        'NS': ('NS', 'Nova Scotia'),
        'NY': ('NY', 'New York'),
        'OH': ('OH', 'Ohio'),
        'OK': ('OK', 'Oklahoma'),
        'ON': ('ON', 'Ontario'),
        'OR': ('OR', 'Oregon'),
        'PA': ('PA', 'Pennsylvania'),
        'PR': ('PR', 'Puerto Rico'),
        'SC': ('SC', 'South Carolina'),
        'SK': ('SK', 'Saskatchewan'),
        'TN': ('TN', 'Tennessee'),
        'TX': ('TX', 'Texas'),
        'UT': ('UT', 'Utah'),
        'VA': ('VA', 'Virginia'),
        'VT': ('VT', 'Vermont'),
        'WA': ('WA', 'Washington'),
        'WI': ('WI', 'Wisconsin'),
        'WV': ('WV', 'West Virginia'),
        'WY': ('WY', 'Wyoming'),
    }
    prefix = site_id[:2].upper()
    return state_map.get(prefix, (prefix, ''))


def get_status(style_url):
    """
    Determine site status from styleUrl.
    n1map = Active (A)
    n0map = Inactive (I) -- has end date
    """
    if 'n1map' in style_url:
        return 'A'
    elif 'n0map' in style_url:
        return 'I'
    return 'A'


def parse_kml(kml_file):
    """Parse KML file and return list of site dictionaries."""
    tree = ET.parse(kml_file)
    root = tree.getroot()

    # handle KML namespace
    ns = {'kml': 'http://earth.google.com/kml/2.1'}

    sites = []

    for placemark in root.findall('.//kml:Placemark', ns):
        # get site ID from name tag
        name_el = placemark.find('kml:name', ns)
        site_id = name_el.text.strip() if name_el is not None else ''

        # get style URL for status
        style_el = placemark.find('kml:styleUrl', ns)
        style_url = style_el.text.strip() if style_el is not None else ''
        status = get_status(style_url)

        # parse HTML description table
        desc_el = placemark.find('kml:description', ns)
        desc_html = desc_el.text if desc_el is not None else ''

        parser = TableParser()
        parser.feed(desc_html)
        table_data = parser.data

        # extract fields from parsed table
        site_name  = table_data.get('Site Name', '')
        lat_str    = table_data.get('Latitude', '')
        lon_str    = table_data.get('Longitude', '')
        elev_str   = table_data.get('Elevation', '')
        start_str  = table_data.get('Start Date', '')
        end_str    = table_data.get('End Date', '')

        # parse dates
        start_date = parse_date(start_str)
        end_date   = parse_date(end_str)

        # parse lat/lon/elevation
        try:
            latitude = float(lat_str) if lat_str else ''
        except ValueError:
            latitude = ''

        try:
            longitude = float(lon_str) if lon_str else ''
        except ValueError:
            longitude = ''

        try:
            elevation = int(float(elev_str)) if elev_str else ''
        except ValueError:
            elevation = ''

        # get state info from site ID prefix
        state_abbrev, state_name = get_state_info(site_id)

        sites.append({
            'network':   'AMoN',
            'siteId':    site_id,
            'siteName':  site_name,
            'status':    status,
            'startDate': start_date,
            'stopDate':  end_date,
            'county':    '',           # not available in KML
            'state':     state_abbrev,
            'latitude':  latitude,
            'longitude': longitude,
            'elevation': elevation,
            'stateName': state_name,
            'siteClass': '',           # not available in KML
        })

    return sites


def write_csv(sites, output_file):
    """Write sites list to CSV file."""
    fieldnames = [
        'network', 'siteId', 'siteName', 'status',
        'startDate', 'stopDate', 'county', 'state',
        'latitude', 'longitude', 'elevation',
        'stateName', 'siteClass'
    ]

    with open(output_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sites)

    print(f"Written {len(sites)} sites to {output_file}")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Convert AMoN KML file to CSV format.'
    )
    parser.add_argument('kml_file',
                        help='input KML file')
    parser.add_argument('csv_file',
                        help='output CSV file')

    args = parser.parse_args()


    kml_file = args.kml_file
    csv_file = args.csv_file

    if len(sys.argv) >= 2:
        kml_file = sys.argv[1]
    if len(sys.argv) >= 3:
        csv_file = sys.argv[2]

    print(f"Parsing {kml_file}...")
    sites = parse_kml(kml_file)
    print(f"Found {len(sites)} sites")

    write_csv(sites, csv_file)

    # print first few rows for verification
    print("\nFirst 5 rows:")
    for s in sites[:5]:
        print(f"  {s['siteId']:6s}  {s['siteName'][:40]:40s}  "
              f"{s['status']}  {s['startDate']}  {s['stopDate']:10s}  "
              f"{s['latitude']:8.4f}  {s['longitude']:10.4f}  "
              f"{str(s['elevation']):6s}  {s['stateName']}")
