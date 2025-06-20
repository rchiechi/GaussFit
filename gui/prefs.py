'''
A preferences frame for settings options, dynamically generated from argparse.
'''
import argparse
import tkinter as tk
from tkinter import ttk, Toplevel, Y, BOTTOM, LEFT, RIGHT, END, BOTH, W
from .colors import GREY
from .tooltip import createToolTip

class PreferencesWindow(Toplevel):
    """
    A modal dialog window for setting advanced preferences.
    This window holds the state of the options and communicates changes
    back to the parent window upon being closed by the save button.
    """
    def __init__(self, parent, opts_parser, opts):
        super().__init__(parent)
        self.title("Preferences")
        
        # This flag will be checked by the parent window to see if we saved.
        self.saved = False
        # This attribute holds the options state. It will be modified by the widgets
        # and then read by the parent window after this dialog is closed.
        # self.result_opts = opts
        
        # Pass the parser, the options state, and the save callback to the frame
        self.prefs = PreferencesFrame(self, opts_parser, opts, self.on_save)
        
        self.transient(parent)
        self.protocol("WM_DELETE_WINDOW", self.on_cancel)
        self.grab_set()

    def on_save(self):
        """Callback for the 'Save' button. Sets the flag and closes the dialog."""
        self.saved = True
        self.destroy()

    def on_cancel(self):
        """Callback for closing the window without saving."""
        self.saved = False
        self.destroy()


class PreferencesFrame(ttk.Frame):
    def __init__(self, master, opts_parser, opts, save_callback): # Note the new 'save_callback'
        super().__init__(master)
        self.opts = opts
        self.opts_parser = opts_parser
        self.master.tk_setPalette(background=GREY, activeBackground=GREY)
        self.pack(fill=BOTH, expand=True, padx=5, pady=5)

        # --- Frame Creation ---
        self.main_content = ttk.Frame(self)
        self.main_content.pack(side=tk.TOP, fill=BOTH, expand=True)

        self.left_pane = ttk.Frame(self.main_content)
        self.left_pane.pack(side=LEFT, fill=Y, padx=5, anchor='n')

        self.right_pane = ttk.Frame(self.main_content)
        self.right_pane.pack(side=RIGHT, fill=Y, padx=5, anchor='n')
        
        self.pane_map = {
            'Input & Output': self.left_pane,
            'Data Processing & Filtering': self.left_pane,
            'Plotting & Visualization': self.right_pane,
            'Statistical & Fitting Parameters': self.right_pane,
            'SLM & Clustering': self.right_pane,
        }

        self.__create_dynamic_options()
        # Pass the save callback to the button creator
        self.__create_buttons(save_callback)

    def __create_buttons(self, save_callback):
        button_frame = ttk.Frame(self)
        button_frame.pack(side=BOTTOM, fill=tk.X, pady=(10, 0))
        # The button now calls the callback provided by the parent Toplevel window
        save_button = ttk.Button(button_frame, text="Save and Close", command=save_callback)
        save_button.pack()

    def __create_dynamic_options(self):
        """
        Dynamically generate options by introspecting the argparse opts_parser.
        """
        for group in self.opts_parser._action_groups:
            parent_pane = self.pane_map.get(group.title)
            if not parent_pane:
                continue

            group_frame = ttk.LabelFrame(parent_pane, text=group.title, padding=10)
            group_frame.pack(fill=BOTH, expand=True, pady=5)

            for action in group._group_actions:
                if action.dest == 'help' or any(opt.startswith('--no') for opt in action.option_strings):
                    continue

                frame = ttk.Frame(group_frame)
                frame.pack(fill='x', pady=4)
                
                # Use a more descriptive label, like the first option string
                label_text = action.option_strings[0].lstrip('-')
                if len(label_text) == 1:
                    label_text = action.option_strings[1].lstrip('-')
                label = ttk.Label(frame, text=f"{label_text}:")
                label.pack(side=LEFT, anchor=W, padx=(0, 5))

                current_value = getattr(self.opts, action.dest, action.default)

                if isinstance(action, (argparse._StoreTrueAction, argparse._StoreFalseAction)):
                    var = tk.BooleanVar(value=current_value)
                    widget = ttk.Checkbutton(frame, variable=var)
                    var.trace_add("write", lambda *_, dest=action.dest, v=var: setattr(self.opts, dest, v.get()))
                else:
                    var = tk.StringVar()
                    widget = ttk.Entry(frame, textvariable=var, width=15)
                    
                    if isinstance(action.nargs, int) and action.nargs > 1:
                        var.set(", ".join(map(str, current_value)))
                        var.trace_add("write", lambda *_, dest=action.dest, v=var, type_func=action.type:
                            self.update_list_opt(dest, v, type_func))
                    else:
                        var.set(str(current_value))
                        var.trace_add("write", lambda *_, dest=action.dest, v=var, type_func=action.type:
                            self.update_scalar_opt(dest, v, type_func))

                widget.pack(side=RIGHT)
                createToolTip(widget, action.help)

    def update_scalar_opt(self, dest, var, type_func):
        """Callback to update a single value in the opts namespace."""
        value_str = var.get()
        try:
            new_value = type_func(value_str) if type_func else value_str
            setattr(self.opts, dest, new_value)
        except (ValueError, TypeError):
            pass

    def update_list_opt(self, dest, var, type_func):
        """Callback to update a list value in the opts namespace."""
        value_str = var.get()
        try:
            new_values = [type_func(v.strip()) for v in value_str.split(',')]
            setattr(self.opts, dest, new_values)
        except (ValueError, TypeError):
            pass
