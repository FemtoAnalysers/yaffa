'''
Application to show the status of the batchfarm.
'''

from curses import wrapper

def center_h(stdscr, line, text):
    '''
    Centers a text in the middle of the screen.
    '''

    _, width = stdscr.getmaxyx()
    x = (width - len(text)) // 2
    stdscr.addstr(line, x, text)

def main(stdscr):
    '''
    Main application.
    '''

    # Clear screen
    stdscr.clear()

    center_h(stdscr, 0, "Batchfarm Monitor")

    stdscr.refresh()
    stdscr.getkey()

wrapper(main)
