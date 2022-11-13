import os
import yagmail
from . import meta


def all_files_from_dir(
    filepath=".",
    extension=".png",
    gmail_username=None,
    separate=False,
    contents=None,
    subject="[dreamcoat] all_files_from_dir",
    to=None,
):
    """Send all files in the filepath with the given extension to the given email.

    Parameters
    ----------
    filepath : str, optional
        File path to send files from, by default ".".
    extension : str, optional
        File extension to send, by default ".png".
    gmail_username : _type_, optional
        Username of the sending gmail account, by default None, in which case this is
        taken from the .dreamcoat .dat files or the user is prompted for input.
    separate : bool, optional
        Whether to send files in separate emails (True) or all together (False), by
        default False.
    contents : str, optional
        Any text to go in the main body of the email, by default None.
    subject : str, optional
        The subject of the email, by default "[dreamcoat] all_files_from_dir".
    to : str, optional
        The email address to send to, by default Nonein which case this is taken from
        the .dreamcoat .dat files or the user is prompted for input.
    """
    if gmail_username is None:
        gmail_username = meta.get_dat_data("gmail_username")
    if to is None:
        to = meta.get_dat_data("to")
    yag = yagmail.SMTP(gmail_username)
    if not filepath.endswith(os.sep):
        filepath += os.sep
    filenames = [filepath + f for f in os.listdir(filepath) if f.endswith(extension)]
    if len(filenames) > 0:
        if separate:
            for f in filenames:
                yag.send(to=to, contents=contents, subject=subject, attachments=f)
        else:
            yag.send(to=to, contents=contents, subject=subject, attachments=filenames)
    else:
        print("No files of given format found in given filepath!")
