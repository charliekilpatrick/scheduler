"""Google Sheets integration for Nickel schedule logging. Requires pip install scheduler[google]."""

import pickle
import os.path
from ._datetime_parse import parse_month_day_year
from .logging import get_logger

logger = get_logger(__name__)


def gsheets_params(current_sheet, template, token, target_dir='', observers_sheet=''):
    """Return dict: CURRENT_SHEET, OBSERVERS_SHEET, target_dir, TEMPLATE, gsheets_token."""
    return {
        'CURRENT_SHEET': current_sheet,
        'OBSERVERS_SHEET': observers_sheet,
        'target_dir': target_dir,
        'TEMPLATE': template,
        'gsheets_token': token,
    }


def initiate_gsheet(token_file):
    """Build Sheets API service from pickle token file."""
    try:
        from googleapiclient.discovery import build
    except ImportError:
        raise ImportError(
            "Google Sheets requires google-api-python-client and google-auth-oauthlib. "
            "Install with: pip install scheduler[google]"
        ) from None
    creds = None
    if os.path.exists(token_file):
        with open(token_file, 'rb') as f:
            creds = pickle.load(f)
    return build('sheets', 'v4', credentials=creds).spreadsheets()


def check_if_tel_on_date(sheet, observers_sheet, date):
    """Return (observer_name, n) if telescope is on date, else (None, 0). date is datetime."""
    response = sheet.values().get(spreadsheetId=observers_sheet, range='Observers').execute()
    data = response['values'][21:]
    year, month, n = '2017', '', 0
    for row in data:
        if not row:
            continue
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
            current_date = parse_month_day_year(month, day, year)
            if current_date.date() == date.date():
                return (row[3], n) if len(row) > 3 else ('Unknown', n)
    return (None, 0)


def need_to_generate_schedule(sheet, current_log, sheet_name):
    """True if sheet_name tab does not exist in current_log."""
    response = sheet.get(spreadsheetId=current_log).execute()
    if response and 'sheets' in response:
        sheet_id = [s['properties']['sheetId'] for s in response['sheets']
                   if s['properties']['title'] == sheet_name]
        return len(sheet_id) == 0
    return None


def copy_log(sheet, template, log_name, current_log, sheet_name, clobber=False):
    """Copy template tab log_name into current_log as sheet_name. Delete existing if clobber."""
    response = sheet.get(spreadsheetId=current_log).execute()
    new_sheet_id = [s['properties']['sheetId'] for s in response['sheets']
                    if s['properties']['title'] == sheet_name]
    if len(new_sheet_id) == 1:
        if not clobber:
            logger.warning("Sheet=%s already exists.", sheet_name)
            return False
        sheet.batchUpdate(spreadsheetId=current_log,
                          body={'requests': [{'deleteSheet': {'sheetId': new_sheet_id[0]}}]}).execute()
    response = sheet.get(spreadsheetId=template).execute()
    if not response or 'sheets' not in response:
        return False
    sheet_id = [s['properties']['sheetId'] for s in response['sheets'] if s['properties']['title'] == log_name]
    if len(sheet_id) != 1:
        return False
    r = sheet.sheets().copyTo(spreadsheetId=template, sheetId=sheet_id[0],
                              body={'destinationSpreadsheetId': current_log}).execute()
    new_sheet_id = r['sheetId']
    sheet.batchUpdate(spreadsheetId=current_log, body={'requests': [{
        'updateSheetProperties': {'properties': {'title': sheet_name, 'sheetId': new_sheet_id}, 'fields': 'title'}
    }]}).execute()
    return True


def populate_nickel_log(sheet, current_log, sheet_name, observations, date=None, observer=None, start_number=None):
    """Fill sheet tab with date, observer, start_number and observation rows."""
    data = sheet.values().get(spreadsheetId=current_log, range=sheet_name).execute()
    log_data = data['values']
    if date:
        log_data[0][1] = date
    if observer:
        log_data[0][3] = observer
    if start_number:
        log_data[7][0] = str(start_number) + '-'
    log_data = log_data + observations
    sheet.values().update(spreadsheetId=current_log, range=sheet_name,
                          valueInputOption='USER_ENTERED', body={'majorDimension': 'ROWS', 'values': log_data}).execute()
