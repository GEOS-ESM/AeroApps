import sys
import argparse
import datetime as dt

def parse_args(args=None):

    parser = argparse.ArgumentParser()

  # Required

    required = parser.add_argument_group('required arguments')

    required.add_argument(
        '--config', metavar='CONFIG', required=True,
        help='Name of configuration file'
    )

  # Optional

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

    return p_args

def make_dt(dt_string):

    if len(dt_string) == 15: return dt_string

    if len(dt_string) <= 8:
        dt_string += 'T000000'
    else:
        dt_string += '000000'

    return dt_string[0:14]
