try:
    from colorama import init, Fore, Style
except ImportError:
    import sys
    print('Error importing colorama module (try pip3 install colorama)')
    sys.exit()

_all__ = ["YELLOW", "WHITE", "RED", "TEAL", "GREEN", "BLUE", "RS"]


# Setup colors
init(autoreset=True)
YELLOW = Fore.YELLOW
WHITE = Fore.WHITE
RED = Fore.RED
TEAL = Fore.CYAN
GREEN = Fore.GREEN
BLUE = Fore.BLUE
RS = Style.RESET_ALL
