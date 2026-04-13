import subprocess

def run_bash_command(command):
    try:
        # Run the command in a subprocess, capture output
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        # Check if the command was successful
        if result.returncode == 0:
            return result.stdout.strip()
        else:
            return f"Error: Command '{command}' returned non-zero exit status {result.returncode}"
    except Exception as e:
        return f"Error: {e}"

# Example usage
command_output = run_bash_command("ls -l")
print("Output of 'ls -l' command:")
print(command_output)