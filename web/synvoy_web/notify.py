"""The one e-mail a visitor may ask for: "your search has finished".

Off unless ``smtp_host``, ``smtp_from`` and ``public_base_url`` are set. Sending never
raises: a mail server that is down must not turn a finished search into a failed one.
The address is deleted from the database right after the attempt.
"""

import smtplib
import ssl
from email.message import EmailMessage


def result_url(settings, job_id):
    return "%s/result.html?id=%s" % (settings.public_base_url.rstrip("/"), job_id)


def compose(settings, to, job_id, title, state, error_text=""):
    message = EmailMessage()
    message["From"] = settings.smtp_from
    message["To"] = to
    link = result_url(settings, job_id)
    if state == "done":
        message["Subject"] = "Your SynVoy search has finished"
        lines = ["Your search has finished:", "", "    " + title, "",
                 "Result: " + link, "",
                 "The result is kept for %d days." % settings.retention_days]
    else:
        message["Subject"] = "Your SynVoy search stopped with an error"
        lines = ["Your search could not be completed:", "", "    " + title, ""]
        if error_text:
            lines += [error_text, ""]
        lines += ["Details: " + link]
    lines += ["", "You get this message because this address was entered when the search",
              "was started. It has now been deleted from the server."]
    message.set_content("\n".join(lines) + "\n")
    return message


def send(settings, message, log=None):
    """True when the mail server accepted the message."""
    try:
        if settings.smtp_port == 465:
            server = smtplib.SMTP_SSL(settings.smtp_host, settings.smtp_port, timeout=20,
                                      context=ssl.create_default_context())
        else:
            server = smtplib.SMTP(settings.smtp_host, settings.smtp_port, timeout=20)
        try:
            if settings.smtp_port != 465 and settings.smtp_starttls:
                server.starttls(context=ssl.create_default_context())
            if settings.smtp_user:
                server.login(settings.smtp_user, settings.smtp_password)
            server.send_message(message)
        finally:
            try:
                server.quit()
            except (smtplib.SMTPException, OSError):
                pass
        return True
    except (smtplib.SMTPException, OSError, ValueError) as exc:
        if log:
            log("e-mail not sent: %s" % exc)
        return False
