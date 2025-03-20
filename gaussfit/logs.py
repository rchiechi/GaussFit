import logging
import colorama as cm

cm.init()

class GaussfitFormatter(logging.Formatter):
    #  https://stackoverflow.com/questions/14844970/modifying-logging-message-format-based-on-message-logging-level-in-python3
    err_fmt = f"{cm.Fore.WHITE}{cm.Back.BLACK}%(asctime)s "\
              f"{cm.Style.BRIGHT}{cm.Fore.CYAN}%(name)s{cm.Fore.WHITE} [{cm.Fore.RED}ERROR{cm.Fore.WHITE}] %(message)s{cm.Style.RESET_ALL}"

    dbg_fmt = f"{cm.Fore.WHITE}{cm.Back.BLACK}%(asctime)s "\
              f"{cm.Style.BRIGHT}{cm.Fore.CYAN}%(name)s{cm.Fore.WHITE} [{cm.Fore.BLUE}DEBUG{cm.Fore.WHITE}] {cm.Style.NORMAL}%(message)s{cm.Style.RESET_ALL}"

    info_fmt = f"{cm.Fore.WHITE}{cm.Back.BLACK}%(asctime)s "\
               f"{cm.Style.BRIGHT}{cm.Fore.CYAN}%(name)s{cm.Fore.WHITE} [{cm.Fore.GREEN}INFO{cm.Fore.WHITE}] %(message)s{cm.Style.RESET_ALL}"

    wrn_fmt = f"{cm.Fore.WHITE}{cm.Back.BLACK}%(asctime)s "\
              f"{cm.Style.BRIGHT}{cm.Fore.CYAN}%(name)s{cm.Fore.WHITE} [{cm.Fore.YELLOW}WARNING{cm.Fore.WHITE}] %(message)s{cm.Style.RESET_ALL}"
    # %(asctime)s - %(name)s - [%(levelname)s] %(message)s'

    def __init__(self):
        super().__init__(fmt='%(asctime)s - %(name)s - [%(levelname)s] %(message)s', datefmt='%H:%M:%S', style='%')

    def format(self, record):

        # Save the original format configured by the user
        # when the logger formatter was instantiated
        format_orig = self._style._fmt

        # Replace the original format with one customized by logging level
        if record.levelno == logging.DEBUG:
            self._style._fmt = GaussfitFormatter.dbg_fmt

        elif record.levelno == logging.INFO:
            self._style._fmt = GaussfitFormatter.info_fmt

        elif record.levelno == logging.WARNING:
            self._style._fmt = GaussfitFormatter.wrn_fmt

        elif record.levelno == logging.ERROR:
            self._style._fmt = GaussfitFormatter.err_fmt

        # Call the original formatter class to do the grunt work
        result = logging.Formatter.format(self, record)
        # Restore the original format configured by the user
        self._style._fmt = format_orig

        return result