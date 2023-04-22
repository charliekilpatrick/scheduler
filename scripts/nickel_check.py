from Logs import *
import sys
from dateutil.parser import parse

date = parse(sys.argv[1])
params = gsheets_params('nickel')
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

