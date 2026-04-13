import subprocess
import time
from datetime import datetime
from asciimatics.screen import Screen
from asciimatics.renderers import BarChart
import curses

# Function to get the number of running and pending jobs using squeue
def get_slurm_jobs(user_only=True):
    cmd = ['squeue', '--noheader', '--format=%T', '-u', USERNAME] if user_only else ['squeue', '--noheader', '--format=%T']
    result = subprocess.run(cmd, capture_output=True, text=True)
    jobs = result.stdout.splitlines()
    running_jobs = jobs.count('R')
    pending_jobs = jobs.count('PD')
    return running_jobs, pending_jobs

# Function to update curses UI
def update_ui(stdscr, total_jobs, user_jobs):
    stdscr.clear()
    stdscr.addstr(0, 0, f"Monitoring Slurm Jobs for User: {USERNAME}")
    stdscr.addstr(1, 0, "Press 'q' to exit.")
    
    stdscr.addstr(3, 0, f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    
    # Display total jobs
    stdscr.addstr(5, 0, "Overall Cluster Stats")
    stdscr.addstr(6, 0, f"Running Jobs: {total_jobs[0]}")
    stdscr.addstr(7, 0, f"Pending Jobs: {total_jobs[1]}")
    
    # Display user's jobs
    stdscr.addstr(9, 0, "User Stats")
    stdscr.addstr(10, 0, f"Running Jobs: {user_jobs[0]}")
    stdscr.addstr(11, 0, f"Pending Jobs: {user_jobs[1]}")
    
    stdscr.refresh()

# Function to render the plot inside the terminal using asciimatics
def render_plot(screen, running_jobs_total, running_jobs_user):
    # Define the size of the plot in the terminal
    max_val = max(running_jobs_total + running_jobs_user) if running_jobs_total and running_jobs_user else 1
    labels = ["Total Running Jobs", "User Running Jobs"]
    data = [running_jobs_total[-20:], running_jobs_user[-20:]]

    # Create a simple bar chart in the terminal
    chart = BarChart(
        screen.height - 4,
        screen.width - 20,
        data,
        char="=",
        labels=labels,
        axes=0,
        scale=max_val,
        border=False,
        colour=[Screen.COLOUR_GREEN, Screen.COLOUR_CYAN],
        bg=Screen.COLOUR_BLACK,
    )

    # Clear the screen and display the plot
    screen.clear_buffer(7, 0, 0)
    screen.print_at("SLURM Job Monitoring - Live Plot", 0, 0, Screen.COLOUR_WHITE, Screen.A_BOLD)
    screen.render_to_buffer(chart, 1, 3)
    screen.refresh()

# Main function to handle ncurses and asciimatics
def monitor_jobs(screen):
    # Initialize data collection
    running_jobs_total = []
    running_jobs_user = []
    timestamps = []

    while True:
        # Get jobs for all users
        total_running, total_pending = get_slurm_jobs(user_only=False)
        user_running, user_pending = get_slurm_jobs(user_only=True)

        # Record data for plotting
        running_jobs_total.append(total_running)
        running_jobs_user.append(user_running)
        timestamps.append(datetime.now().strftime('%H:%M:%S'))

        # Update the plot using asciimatics
        render_plot(screen, running_jobs_total, running_jobs_user)

        # Wait for 1 second before refreshing
        time.sleep(1)

# Wrapper function to manage the curses-based UI for job stats
def curses_ui(stdscr):
    USERNAME = subprocess.run(['whoami'], capture_output=True, text=True).stdout.strip()

    while True:
        total_running, total_pending = get_slurm_jobs(user_only=False)
        user_running, user_pending = get_slurm_jobs(user_only=True)
        update_ui(stdscr, (total_running, total_pending), (user_running, user_pending))

        stdscr.timeout(1000)  # Set a timeout to allow the curses UI to refresh
        key = stdscr.getch()
        if key == ord('q'):
            break  # Exit the loop if 'q' is pressed

if __name__ == "__main__":
    USERNAME = subprocess.run(['whoami'], capture_output=True, text=True).stdout.strip()

    # Launch the job monitor in two threads: one for curses, one for plotting with asciimatics
    try:
        # Use two threads or subprocesses if needed; however, for now, we will alternate UI tasks
        curses.wrapper(curses_ui)  # Manages job stats display
        Screen.wrapper(monitor_jobs)  # Manages the plot display
    except KeyboardInterrupt:
        pass  # Handle Ctrl+C gracefully

