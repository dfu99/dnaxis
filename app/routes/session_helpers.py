from functools import wraps
from flask import session, redirect, url_for, flash


def require_session_keys(*keys):
    """Decorator that redirects to home if any required session keys are missing."""
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            for key in keys:
                if key not in session:
                    flash('Please start from the beginning.')
                    return redirect(url_for('main.home'))
            return f(*args, **kwargs)
        return decorated_function
    return decorator
