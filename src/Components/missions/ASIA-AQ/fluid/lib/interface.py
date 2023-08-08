
import os
import sys
import argparse
import datetime as dt

def parse_args(args=None):

    parser = argparse.ArgumentParser()

#   # required
#   required = parser.add_argument_group('required arguments')

#   required.add_argument(
#       '--rc', metavar='RESOURCE FILE', required=True,
#       help='Name of resource file'
#   )

    # optional
    parser.add_argument(
        '--rc', metavar='RESOURCE', default='',
        help='Name of resource file (default: %(default)s)'
    )
    parser.add_argument(
        '--config', metavar='CONFIG', default=[],action='append',
        help='Name of configuration file or directory'
    )
    parser.add_argument(
        '--reset', metavar='RESET', default=[],action='append',
        help='Name of resources to reset'
    )
    parser.add_argument(
        '--motif', metavar='MOTIF', default=[],action='append',
        help='Name of motif file or directory'
    )
    parser.add_argument(
        '--theme', metavar='THEME', default=[],action='append',
        help='Name of configuration file or directory referencing a theme'
    )
    parser.add_argument(
        '-g', '--geometry', metavar='GEOMETRY', default='1024x768',
        help='Image size in pixels (default: %(default)s)'
    )
    parser.add_argument(
        '-r', '--region', metavar='REGION', default='',
        help='Name of region (default: %(default)s)'
    )
    parser.add_argument(
        '-f', '--field', metavar='FIELD', default='',
        help='Name of field (default: %(default)s)'
    )
    parser.add_argument(
        '-p', '--plot', metavar='PLOT', default='',
        help='Name of field (default: %(default)s)'
    )
    parser.add_argument(
        '-l', '--level', metavar='LEVEL', default='',
        help='Pressure level (default: %(default)s)'
    )
    parser.add_argument(
        '-s', '--stream', metavar='STREAM', default='G5FPFC',
        help='Name of stream (default: %(default)s)'
    )
    parser.add_argument(
        '-c', '--collection', metavar='COLLECTION', default='',
        help='Name of collection (default: %(default)s)'
    )
    parser.add_argument(
        '--layer', metavar='LAYER', default='',
        help='Name of layer (default: %(default)s)'
    )
    parser.add_argument(
        '-t', '--time_dt', metavar='YYYYMMDDTHHMMSS', default=None,
        help='Time in ISO format'
    )
    parser.add_argument(
        '--fcst_dt', metavar='YYYYMMDDTHHMMSS', default=None,
        help='Forecast start time in ISO format'
    )
    parser.add_argument(
        '--start_dt', metavar='YYYYMMDDTHHMMSS', default=None,
        help='Start time in ISO format'
    )
    parser.add_argument(
        '--end_dt', metavar='YYYYMMDDTHHMMSS', default=None,
        help='Ending time in ISO format'
    )
    parser.add_argument(
        '--t_deltat', metavar='HOURS', type=int, default='3',
        help='Time increment in hours (default: %(default)s)'
    )
    parser.add_argument(
        '-o', '--oname', metavar='ONAME', default='%Y%m%dT%H%M%S.png',
        help='Output filename (default: %(default)s)'
    )
    parser.add_argument(
        '-m', '--match', metavar='MATCH', default='%d:01,06,11,16,21,26',
        help='Output filename (default: %(default)s)'
    )
    parser.add_argument(
        '-b', '--basemap', metavar='BASEMAP', default=None,
        help='Basemap name (default: %(default)s)'
    )

    parser.add_argument(
        '--no_label', action='store_true', help='Label on/off flag'
    )

    parser.add_argument(
        '--no_logo', action='store_true', help='Logo on/off flag'
    )

    parser.add_argument(
        '--no_title', action='store_true', help='Title on/off flag'
    )

    parser.add_argument(
        '--plot_only', action='store_true', help='Plot only flag'
    )

    parser.add_argument(
        '--lights_off', action='store_true', help='Lights on/off flag'
    )

    parser.add_argument(
        '--label_size', metavar='LABEL_SIZE', default=None,
        help='Size of lat/lon labels (default: %(default)s)'
    )

    parser.add_argument(
        '--tick_label_size', metavar='TICK_LABEL_SIZE', default=None,
        help='Size of tick labels (default: %(default)s)'
    )

    parser.add_argument(
        '--fields', metavar='FIELDS', default='',
        help='Comma-separated list of fields (default: %(default)s)'
    )

    parser.add_argument(
        '--exclude', metavar='FIELDS', default='',
        help='Comma-separated list of fields (default: %(default)s)'
    )

    parser.add_argument(
        '--navigate', metavar='NAVIGATE', default='on',
        help='navigate on/off (default: %(default)s)'
    )


    if not len(args):
        parser.print_help()
        sys.exit(1)

    p_args = vars(parser.parse_args(args))

    if p_args['time_dt'] is None:
        now    = dt.datetime.utcnow()
        hour   = int(now.hour / 12) * 12
        hhmmss = "%06d"%(hour*10000,)
        p_args['time_dt'] = now.strftime('%Y%m%dT' + hhmmss)

    p_args['time_dt'] = make_dt(p_args['time_dt'])

    if p_args['start_dt'] is None:
        p_args['start_dt'] = make_dt(p_args['time_dt'])
    else:
        p_args['start_dt']  = make_dt(p_args['start_dt'])

    if p_args['end_dt'] is None:
        p_args['end_dt'] = make_dt(p_args['start_dt'])
    else:
        p_args['end_dt']  = make_dt(p_args['end_dt'])

    if p_args['fcst_dt']:
        p_args['fcst_dt'] = make_dt(p_args['fcst_dt'])
        p_args['fcst_dt'] = dt.datetime.strptime(p_args['fcst_dt'],'%Y%m%dT%H%M%S')

    p_args['time_dt']  = dt.datetime.strptime(p_args['time_dt'],'%Y%m%dT%H%M%S')
    p_args['start_dt'] = dt.datetime.strptime(p_args['start_dt'],'%Y%m%dT%H%M%S')
    p_args['end_dt']   = dt.datetime.strptime(p_args['end_dt'],'%Y%m%dT%H%M%S')
    p_args['t_deltat'] = dt.timedelta(hours=p_args['t_deltat'])

#   Record user activity

    s_args = [ v.split('-')[-1] for v in args if str(v)[0] == '-' ]
    s_args = set(s_args)

    user = {}
    if {'region', 'r'} & s_args: user['region'] = 1
    if {'field',  'f'} & s_args: user['field'] = 1
    if {'plot',   'p'} & s_args: user['plot'] = 1
    if {'level',  'l'} & s_args: user['level'] = 1
    if {'stream', 's'} & s_args: user['stream'] = 1
    if {'oname',  'o'} & s_args: user['oname'] = 1
    if {'collection', 'c'} & s_args: user['collection'] = 1
    if {'time_dt'} & s_args: user['time_dt'] = 1
    if {'start_dt'} & s_args: user['start_dt'] = 1
    if {'end_dt'} & s_args: user['end_dt'] = 1
    if {'t_deltat'} & s_args: user['t_deltat'] = 1

    if user.get('plot', 0):
        user['field'] = 1
        p_args['field'] = p_args['plot']

    p_args['user'] = user

    return p_args

def make_dt(dt_string):

    if len(dt_string) == 15: return dt_string

    if len(dt_string) <= 8:
        dt_string += 'T000000'
    else:
        dt_string += '000000'

    return dt_string[0:14]
