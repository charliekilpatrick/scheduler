from scheduler.Logs import *
import sys
import os
from dateutil.parser import parse

current_sheet = os.environ['NICKEL_SCHEDULE_SHEET']
observers_sheet = os.environ['NICKEL_OBSERVERS_SHEET']
template_sheet = os.environ['NICKEL_TEMPLATE_SHEET']
gsheet_token = os.environ['GSHEETS_TOKEN']

date = parse(sys.argv[1])
params = gsheets_params(current_sheet, template_sheet, gsheet_token,
    observers_sheet=observers_sheet)
sheet = initiate_gsheet(params['gsheets_token'])
observers_sheet = params['OBSERVERS_SHEET']
todays_sheet_name = date.strftime('%Y%m%d')

(observer, n) = check_if_tel_on_date(sheet, observers_sheet, date)

if not observer:
    print(0)
else:
    check = need_to_generate_schedule(sheet, params['CURRENT_SHEET'],
        todays_sheet_name)
    if check:
        print(1)
    else:
        print(0)

