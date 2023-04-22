import pickle
import os.path, sys
from googleapiclient.discovery import build
from google_auth_oauthlib.flow import InstalledAppFlow
from google.auth.transport.requests import Request
from dateutil.parser import parse

def gsheets_params(current_sheet, template, token, target_dir='', 
    observers_sheet=''):

    params = {
        # This is the google sheet that holds all of the logs
        'CURRENT_SHEET':current_sheet,

        # This is the google sheet for recording observers on a particular date
        'OBSERVERS_SHEET':observers_sheet,

        # This is the local directory where target lists will be stored
        'target_dir': target_dir, 

        # Now input general params that are same for all telescopes
        # This is the template google sheet that holds blank copies of logs
        'TEMPLATE': template,

        # This is the token for needed for accessing google sheets API
        'gsheets_token': token,
    }

    return(params)

def initiate_gsheet(token_file):
    creds = None
    if os.path.exists(token_file):
        with open(token_file, 'rb') as token:
                creds = pickle.load(token)

    service = build('sheets', 'v4', credentials=creds)
    sheet = service.spreadsheets()

    return(sheet)

# Check if input date is an observing date for tel and if so return observer
# This assumes standard formatting for observers google sheet (see sheets above)
# Date should be a datetime object with the date we want to check
def check_if_tel_on_date(sheet, observers_sheet, date):
    response = sheet.values().get(spreadsheetId=observers_sheet,
        range='Observers').execute()

    data = response['values'][21:]

    year = '2017'    # Year in the sheet
    month = ''   # Month in the sheet
    day = ''     # Day in the sheet
    n = 0        # This day is the nth day of the quarter
    for row in data:
        if row:
            date_value = row[0]
            if 'Q' in date_value and date_value.startswith('20'):
                n = 0
                year = date_value.split('Q')[0]
            elif date_value.isdigit() and int(date_value) > 100:
                year = date_value
            elif not date_value.isdigit() and not date_value.startswith('20'):
                month = date_value
                continue
            else:
                n += 1
                day = date_value
                current_date = parse(month + ' ' + day + ', ' + year)

                if current_date.date() == date.date():
                    if len(row)>3:
                        return(row[3], n)
                    else:
                        return('Unknown', n)

    return(None, 0)


def need_to_generate_schedule(sheet, current_log, sheet_name):
    response = sheet.get(spreadsheetId=current_log).execute()

    if response and 'sheets' in response.keys():
        sheet_id = [sheet['properties']['sheetId']
            for sheet in response['sheets']
            if sheet['properties']['title'] == sheet_name]

        if len(sheet_id)==0:
            return(True)
        else:
            return(False)

    else:
        return(None)

# Copies the blank template log into the input google sheet.  Need to have the
# log_name from the correct tab on the google sheet with templates.  The tab
# will be renamed to sheet_name in the google sheet current_log.
def copy_log(sheet, template, log_name, current_log, sheet_name, clobber=False):
    # First do a check if there exists a sheet with name "sheet_name"
    # in current sheet.  If yes and clobber=False, then return(False).
    # If yes and clobber=True, proceed and re-copy blank sheet into
    # log then return(True).  If no, proceed normally.

    response = sheet.get(spreadsheetId=current_log).execute()

    new_sheet_id = [sheet['properties']['sheetId']
            for sheet in response['sheets']
            if sheet['properties']['title'] == sheet_name]


    if len(new_sheet_id)==1:
        if not clobber:
            warning = 'WARNING: sheet={s} already exists!\n'
            warning += 'Exiting...'
            print(warning.format(s=sheet_name))
            return(False)
        else:
            warning = 'WARNING: sheet={s} already exists and clobber=True\n'
            warning += 'Clobbering...\n\n'
            print(warning.format(s=sheet_name))

            # Since we're clobbering, delete the old sheet and start over
            body = {'requests': [{'deleteSheet': {'sheetId': new_sheet_id[0]}}]}
            req = sheet.batchUpdate(spreadsheetId=current_log, body=body)
            response = req.execute()
    else:
        new_sheet_id = None

    response = sheet.get(spreadsheetId=template).execute()

    if response and 'sheets' in response.keys():
        sheet_id = [sheet['properties']['sheetId']
            for sheet in response['sheets']
            if sheet['properties']['title'] == log_name]

        if len(sheet_id)==1:
            sheet_id = sheet_id[0]
            body={'destinationSpreadsheetId': current_log}

            req = sheet.sheets().copyTo(spreadsheetId=template,
                sheetId=sheet_id, body=body)
            response = req.execute()

            new_sheet_id = response['sheetId']
            new_name = 'Copy of '+log_name

            body = {'requests': [{
                'updateSheetProperties': {
                    'properties': {
                        'title': sheet_name, 'sheetId': new_sheet_id
                    },
                    'fields': 'title'
                }
            }]}

            req = sheet.batchUpdate(spreadsheetId=current_log, body=body)
            response = req.execute()

            return(True)

        else:
            return(False)
    else:
        return(False)

    return(False)

def populate_nickel_log(sheet, current_log, sheet_name, observations,
    date=None, observer=None, start_number=None):

    # Grab data from the current empty sheet
    data = sheet.values().get(spreadsheetId=current_log,
        range=sheet_name).execute()

    log_data = data['values']

    # First populate date-specific metadata
    if date:
        log_data[0][1] = date
    if observer:
        log_data[0][3] = observer
    if start_number:
        log_data[7][0] = str(start_number) + '-'

    # Append the list of observations to the log
    log_data = log_data + observations

    body = {'majorDimension': 'ROWS', 'values': log_data}

    result = sheet.values().update(spreadsheetId=current_log, range=sheet_name,
        valueInputOption='USER_ENTERED', body=body).execute()


