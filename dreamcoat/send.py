import os
import yagmail
from . import meta


def one_file(
    filename,
    filepath=".",
    gmail_username=None,
    gmail_password=None,
    contents=None,
    subject=None,
    to=None,
):
    """Send all files in the filepath with the given extension to the given email.

    Parameters
    ----------
    filename : str
        The name of the file to send.
    filepath : str, optional
        File path to send files from, by default ".".
    gmail_username : str, optional
        Username of the sending gmail account, by default None, in which case this is
        taken from the .dreamcoat .dat files or the user is prompted for input.
    gmail_password : str, optional
        Password of the sending gmail account, by default None, in which case it is
        requested from the keyring or asked of the user.
    contents : str, optional
        Any text to go in the main body of the email, by default None.
    subject : str, optional
        The subject of the email, by default "[dreamcoat] {filepath}{filename}".
    to : str, optional
        The email address to send to, by default Nonein which case this is taken from
        the .dreamcoat .dat files or the user is prompted for input.
    """
    if gmail_username is None:
        gmail_username = meta.get_dat_data("gmail_username")
    if gmail_password is None:
        gmail_password = meta.get_dat_data("gmail_password")
    if to is None:
        to = meta.get_dat_data("to")
    if gmail_password is None:
        yag = yagmail.SMTP(gmail_username)
    else:
        yag = yagmail.SMTP(gmail_username, gmail_password)
    if not filepath.endswith(os.sep):
        filepath += os.sep
    if subject is None:
        subject = "[dreamcoat] {}{}".format(filepath, filename)
    f = filepath + filename
    yag.send(to=to, contents=contents, subject=subject, attachments=f)


def all_files_from_dir(
    filepath=".",
    extension="zip",
    gmail_username=None,
    gmail_password=None,
    separate=False,
    contents=None,
    subject=None,
    to=None,
):
    """Send all files in the filepath with the given extension to the given email.

    Parameters
    ----------
    filepath : str, optional
        File path to send files from, by default ".".
    extension : str, optional
        File extension to send, by default "zip".
    gmail_username : str, optional
        Username of the sending gmail account, by default None, in which case this is
        taken from the .dreamcoat .dat files or the user is prompted for input.
    gmail_password : str, optional
        Password of the sending gmail account, by default None, in which case it is
        requested from the keyring or asked of the user.
    separate : bool, optional
        Whether to send files in separate emails (True) or all together (False), by
        default False.
    contents : str, optional
        Any text to go in the main body of the email, by default None.
    subject : str, optional
        The subject of the email, by default "[dreamcoat] {filepath}*.{extension}".
    to : str, optional
        The email address to send to, by default Nonein which case this is taken from
        the .dreamcoat .dat files or the user is prompted for input.
    """
    if gmail_username is None:
        gmail_username = meta.get_dat_data("gmail_username")
    if gmail_password is None:
        gmail_password = meta.get_dat_data("gmail_password")
    if to is None:
        to = meta.get_dat_data("to")
    if gmail_password is None:
        yag = yagmail.SMTP(gmail_username)
    else:
        yag = yagmail.SMTP(gmail_username, gmail_password)
    if not filepath.endswith(os.sep):
        filepath += os.sep
    filenames = [
        filepath + f for f in os.listdir(filepath) if f.endswith("." + extension)
    ]
    if subject is None:
        subject = "[dreamcoat] {}*.{}".format(filepath, extension)
    if len(filenames) > 0:
        if separate:
            for f in filenames:
                yag.send(to=to, contents=contents, subject=subject, attachments=f)
        else:
            yag.send(to=to, contents=contents, subject=subject, attachments=filenames)
    else:
        print("No files of given format found in given filepath!")


def delete_from_dir(filepath=".", extension="zip"):
    """Delete all files with a given extension from a given directory.

    Parameters
    ----------
    filepath : str, optional
        File path to send files from, by default ".".
    extension : str, optional
        File extension to send, by default "zip".
    """
    if not filepath.endswith(os.sep):
        filepath += os.sep
    filenames = [
        filepath + f for f in os.listdir(filepath) if f.endswith("." + extension)
    ]
    if len(filenames) > 0:
        for f in filenames:
            print(f)
            os.remove(f)
