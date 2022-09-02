# TODO
# Change into a class?
# Change import from own sources

# Rewrite it so that I use a authenticatet ldap user and tries the user name and password from that user and see I can create A
# bind

from functools import wraps
from flask import session, flash, redirect

# ldap stuff
from ldap3 import Server, Connection, ALL, SUBTREE
from ldap3.core.exceptions import LDAPException, LDAPBindError

# import own stuff
# from ../config import LDAP_SERVER, LDAP_BASE_DN


# works as a decorator
def logged_in_required(func):
    @wraps(func)
    def inner(*args, **kwargs):
        if "active_user" in session:
            return func(*args, **kwargs)
        else:
            flash("Please log in first", "danger")
            return redirect("login")

    return inner


def ldap_authentication(user_name: str, user_pwd: str) -> bool:
    """
    Modified from https://soshace.com/integrate-ldap-authentication-with-flask/
    """
    ldap_user_name = user_name.strip()
    ldap_user_pwd = user_pwd.strip()
    user = f"cn={ldap_user_name},{LDAP_BASE_DN}"
    server = Server(LDAP_SERVER, get_info=ALL)
    connection = Connection(server, user=user, password=ldap_user_pwd)

    return connection.bind()
