from flask import Flask, render_template, url_for, request, flash, redirect, session

# init app and objects
app = Flask(__name__)

# config
app.secret_key = '70c6a968ed1ada341dbcbf252b3ea3cf'

# ldap-login 
@app.before_request
def make_session_expire():
    '''
    Makes the session expire after 5 minutes, if no GET REQUEST is made.
    '''
    session.permanent = True
    app.permanent_session_lifetime = timedelta(minutes=5)
    
