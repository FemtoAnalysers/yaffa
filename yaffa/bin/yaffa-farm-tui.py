'''
Application to show the status of the batchfarm.
'''

import time
import subprocess
import curses
from curses import wrapper

def center_h(stdscr, line, text):
    '''
    Centers a text in the middle of the screen.
    '''

    _, width = stdscr.getmaxyx()
    x = (width - len(text)) // 2
    stdscr.addstr(line, x, text)

def output_of_cmd(command):
    '''
    Returns the output of a shell command
    '''

    return subprocess.run(command, shell=True, check=True, stdout=subprocess.PIPE).stdout.decode().strip()

def main(stdscr):
    '''
    Main application.
    '''
    curses.curs_set(0)

    lines, cols = stdscr.getmaxyx()
    while True:
        # Check if the dimensions have changed, if so clear the screen
        if (lines, cols) != stdscr.getmaxyx():
            lines, cols = stdscr.getmaxyx()
            stdscr.clear()

        center_h(stdscr, 0, "Batchfarm Monitor")

        winInfo = curses.newwin(6, 162 // 2, 1, 0)
        winInfo.border()
        winInfo.addstr("Status ", curses.A_BOLD)

        # All jobs
        totalJobs = output_of_cmd('squeue -h | wc -l')
        totalRunningJobs = output_of_cmd('squeue -h -t R | wc -l')
        totalPendingJobs = output_of_cmd('squeue -h -t PD | wc -l')
        totalHeldJobs = output_of_cmd('squeue -h | grep JobHeldUser | wc -l')

        winInfo.addstr(1, 2, f'Total jobs:   {totalJobs:>6} │ ')
        winInfo.addstr(2, 2, f'Running jobs: {totalRunningJobs:>6} │ ')
        winInfo.addstr(3, 2, f'Pending jobs: {totalPendingJobs:>6} │ ')
        winInfo.addstr(4, 2, f'Held jobs:    {totalHeldJobs:>6} │ ')

        # My jobs
        myJobs = output_of_cmd('squeue --me -h | wc -l')
        myRunningJobs = output_of_cmd('squeue --me -h -t R | wc -l')
        myPendingJobs = output_of_cmd('squeue --me -h -t PD | wc -l')
        myHeldJobs = output_of_cmd('squeue --me -h | grep JobHeldUser | wc -l')

        winInfo.addstr(1, 25, f'Your jobs:         {myJobs:>6} │ ')
        winInfo.addstr(2, 25, f'Your running jobs: {myRunningJobs:>6} │ ')
        winInfo.addstr(3, 25, f'Your pending jobs: {myPendingJobs:>6} │ ')
        winInfo.addstr(4, 25, f'Your held jobs:    {myHeldJobs:>6} │ ')

        # Their jobs
        theirJobs = output_of_cmd('squeue -h | grep -v $(whoami) | wc -l')
        theirRunningJobs = output_of_cmd('squeue -h -t R | grep -v $(whoami) | wc -l')
        theirPendingJobs = output_of_cmd('squeue -h -t PD | grep -v $(whoami) | wc -l')
        theirHeldJobs = output_of_cmd('squeue -h | grep -v $(whoami) | grep JobHeldUser | wc -l')

        winInfo.addstr(1, 53, f'Their jobs:         {theirJobs:>6}')
        winInfo.addstr(2, 53, f'Their running jobs: {theirRunningJobs:>6}')
        winInfo.addstr(3, 53, f'Their pending jobs: {theirPendingJobs:>6}')
        winInfo.addstr(4, 53, f'Their held jobs:    {theirHeldJobs:>6}')

        winInfo.refresh()

        stdscr.refresh()
        time.sleep(1)

wrapper(main)
