#!/usr/bin/python3

__all__ = ["YELLOW", "WHITE", "RED", "TEAL", "GREEN", "BLUE", "RS"]

from colorama import init,Fore,Back,Style

# Setup colors
init(autoreset=True)
YELLOW=Fore.YELLOW
WHITE=Fore.WHITE
RED=Fore.WHITE
TEAL=Fore.CYAN
GREEN=Fore.GREEN
BLUE=Fore.BLUE
RS=Style.RESET_ALL

#YELLOW="\033[1;33m"
#WHITE="\033[0m"
#RED="\033[1;31m"
#TEAL="\033[1;36m"
#GREEN="\033[1;32m"
#BLUE="\033[1;34m"
#RS="\033[0m"
#CL="\033[2K"
