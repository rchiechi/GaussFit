import tkinter as tk
from tkinter import ttk, Toplevel, W, E, N, S, BOTH, LEFT, RIGHT, VERTICAL, Y

class AdvancedSettingsWindow(Toplevel):
    def __init__(self, master, opts, parser):
        super().__init__(master)
        self.title("Advanced Settings")
        self.opts = opts
        self.parser = parser

        self.transient(master)
        self.grab_set()

        self.main_frame = ttk.Frame(self, padding="10")
        self.main_frame.pack(expand=True, fill=BOTH)

        self.create_widgets()
        self.protocol("WM_DELETE_WINDOW", self.destroy)

    def create_widgets(self):
        # We can group arguments for clarity
        # The parser object from argparse is available at self.parser
        for group in self.parser._action_groups:
            if group.title in ['positional arguments', 'optional arguments']:
                continue # Skip default groups if they are not needed

            group_frame = ttk.LabelFrame(self.main_frame, text=group.title, padding="10")
            group_frame.pack(fill=BOTH, expand=True, padx=5, pady=5)

            for action in group._group_actions:
                # Skip help action
                if action.dest == 'help':
                    continue

                frame = ttk.Frame(group_frame)
                frame.pack(fill='x', pady=2)

                # Get the current value from the options namespace
                current_value = getattr(self.opts, action.dest, action.default)

                if isinstance(action, argparse._StoreTrueAction):
                    var = tk.BooleanVar(value=current_value)
                    widget = ttk.Checkbutton(frame, text=action.help, variable=var)
                    widget.pack(side=LEFT)
                    # When the checkbox is toggled, update the opts namespace
                    widget.config(command=lambda dest=action.dest, v=var: setattr(self.opts, dest, v.get()))

                elif isinstance(action, argparse._StoreAction):
                    label = ttk.Label(frame, text=f"{action.dest}:", wraplength=200)
                    label.pack(side=LEFT, padx=(0, 5))
                    
                    var = tk.StringVar(value=str(current_value))
                    widget = ttk.Entry(frame, textvariable=var)
                    widget.pack(side=LEFT, fill='x', expand=True)
                    # Update the opts namespace when the entry changes
                    var.trace_add("write", lambda *_, dest=action.dest, v=var: self.update_opt(dest, v, action.type))
                
                # Add a tooltip with the help string
                if action.help:
                     # Using the existing createToolTip function
                     createToolTip(widget, action.help)


    def update_opt(self, dest, var, type_func):
        """Update the options namespace with type conversion."""
        value = var.get()
        try:
            # Convert the value back to its original type (e.g., int, float)
            if type_func:
                setattr(self.opts, dest, type_func(value))
            else:
                setattr(self.opts, dest, value)
        except (ValueError, TypeError):
            # Handle cases where the input isn't the correct type, e.g., show an error.
            # For simplicity, we'll just ignore invalid values here.
            pass