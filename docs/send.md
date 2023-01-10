# Sending data by email

## Sending data

To send data by email with dreamcoat you will first need to set up a Gmail account with an app-specific password.  This account is used to send the emails and its username goes into the `gmail_username` kwarg below.

### One specific file

```python
import dreamcoat as dc

dc.send.one_file(
    filename,
    filepath=".",
    gmail_username=None,
    separate=False,
    contents=None,
    subject="[dreamcoat] all_files_from_dir",
    to=None,
)
```

  * [Yagmail](https://github.com/kootenpv/yagmail) does the work behind the scenes, so see its documentation for more info about the other kwargs.
  * Your `gmail_username` (and password) can optionally be put into the .dreamcoat directory, as can the `to` value.

### All files with given extension

To attach all files with a given path and extension to an email and send it, use

```python
dc.send.all_files_from_dir(
    filepath=".",
    extension="zip",
    gmail_username=None,
    separate=False,
    contents=None,
    subject="[dreamcoat] all_files_from_dir",
    to=None,
)
```

  * Only files matching the specified `extension` are sent.
  * `separate` determines whether to send each file in a separate email or attach them all to the same one.

## Tidy up after sending

After sending the files, you can delete all the files with a given extension from a directory using `delete_from_dir`:

```python
dc.send.delete_from_dir(filepath=".", extension="zip")
```

Note that this will indiscriminately delete all files in the specified `filepath` with the specified `extension` - not just the ones you sent!  Use `os.remove()` if you need to be more selective.
