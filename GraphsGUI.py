import csv
import os
import tkinter as tk
from tkinter import ttk
import matplotlib
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import DataHandler

LARGE_FONT = ("Verdana", 12)


class GraphGUI(tk.Tk):
    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)
        tk.Tk.iconbitmap(self, default='./res/m2micon.ico')
        tk.Tk.wm_title(self, "Moment to moment graph generator")
        container = tk.Frame()
        container.pack(side='top', fill='both', expand=True)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)
        self.handler = DataHandler.DataHandler()
        self.setup_graph()
        self.frames = {}
        self.user_id = self.handler.user_ids[0]
        self.objective_id = self.handler.learn_obj_ids[0]
        for F in (StartPage, GraphPage):
            new_frame = F(container, self)
            self.frames[F] = new_frame
            new_frame.grid(row=0, column=0, sticky='nsew')
        self.show_frame(GraphPage)
        self.method = 'all'

    def show_frame(self, cont):
        frame = self.frames[cont]
        frame.tkraise()

    def update_graph(self, huh=None):
        # cmap = plt.get_cmap('cool')
        graph_n, graph_l, graph_f, o_graph, split, answers = \
            self.handler.get_graph_variables(
            self.user_id, method=self.method, oid=self.objective_id)
        t_list = self.handler.get_graph_variables(self.user_id,
                                                        method=self.method,
                                                        oid=self.objective_id)
        self.a.clear()
        self.a.plot(range(len(graph_n)), graph_n, label="P(Jn)")
        # self.a.plot(range(len(graph_n)), graph_l, label="P(Jl)")
        self.a.plot(range(len(graph_n)), graph_f, label="P(Jf)")
        self.a.plot(range(len(answers)), answers, label="Answers")
        height = max(max(graph_n), max(graph_l), max(graph_f))
        low = min(min(graph_n), min(graph_l), min(graph_f))
        self.a.legend()
        self.a.plot([split[0], split[0]], [low, height], color="black")
        self.a.text(split[0], height, str(split[1]),
                    horizontalalignment='center',
                    verticalalignment='center', bbox=dict(facecolor='white',
                                                          edgecolor='white',
                                                          alpha=1.0))
        boundary_list = self.handler.boundary_list[:]
        color_list = self.handler.color_list[:len(boundary_list) - 1]
        for b1, b2, c in zip(boundary_list[:-1], boundary_list[1:],
                             color_list):
            self.a.broken_barh([(b1, b2 - b1)],
                               (low + .15 * (height-low), height - .15 * (
                                       height-low)),
                               facecolors=c)
        print(o_graph)

    def setup_graph(self):
        self.f = Figure(figsize=(5, 5), dpi=100)
        self.a = self.f.add_subplot(111)
        self.a.plot([1, 1], [1, 0])

    def upd_user(self, user_id):
        self.user_id = user_id
        self.update_graph(0)

    def upd_learn_goal(self, l_o_id):
        self.objective_id = l_o_id
        self.update_graph(0)

    def move_user(self, incr):
        idx = self.handler.get_users().tolist().index(self.user_id) + incr
        self.user_id = self.handler.get_users()[
            idx % len(self.handler.get_users())]
        return self.user_id

    def save_graph(self, fname=None):
        if fname is None:
            fname = 'graphs/student {} objective {}.png'.format(self.user_id,
                                                                self.objective_id)
        self.f.savefig(fname=fname)

    def save_all_graphs(self, dirname='graphs_simone/'):
        if not os.path.isdir(dirname):
            os.makedirs(dirname)
        with open(dirname+'spikes.csv', 'w', newline='') as csv_file:
            writer = csv.writer(csv_file)
            cat = ["Spikiness", "Peaks", "Transitional peaks"]
            kind = ["in general", "voormeting", "instructie",
                    "non-adaptief", "adaptief na les",
                    "adaptief herhaling", "nameting"]
            header_row = ["Student", "Leerdoel"]
            for c in cat:
                for k in kind:
                    header_row.append(c+" "+k)
            writer.writerow(header_row)
            errors = 0
            for user in self.handler.get_users():
                for learn_obj in np.unique(self.handler.learn_obj_ids):
                    try:
                        f = matplotlib.pyplot.figure(figsize=(5, 5), dpi=100)
                        axes = matplotlib.pyplot.gca()
                        axes.set_ylim([-1.15, 1.1])
                        a = f.add_subplot(111)
                        graph_n, graph_l, graph_f, o_graph, \
                            split, answers = \
                            self.handler.get_graph_variables(user,
                                                             method=self.method,
                                                             oid=learn_obj)
                        x = [1+i for i in range(len(graph_n))]
                        # a.plot(range(len(graph_n)), graph_f, label="P(Jf)")
                        height = 1.  # max([max(graph_n), max(graph_l)])
                        low = -1.  # min([min(graph_n), min(graph_l)])
                        height = height + .05*(height-low)
                        a.plot([split[0], split[0]], [low, height], color="black")
                        if split[0] is not None:
                            a.text(max(split[0], len(graph_n)/50), height,
                                   str(split[1]),
                                   horizontalalignment='center',
                                   verticalalignment='center',
                                   bbox=dict(facecolor='white', edgecolor='white',
                                             alpha=1.0))
                        # a.plot(x, graph_n,
                        #        label="Nieuwe curve", color="black")
                        # a.plot(x, graph_l, label="Oude curve",
                        #        color="deeppink", dashes=[2, 2])
                        a.plot(x, o_graph, label="Oude oude curve",
                               color="cyan")
                        new_low = -1.1  # low - .1 * (height - low)
                        new_height = -1.05  # new_low + .05 * (height - low)
                        a.plot(range(1, len(answers)+1),
                               [new_height if a == 1 else new_low for a in
                                answers],
                               color="red", label="Answers")
                        # a.legend()

                        boundary_list = self.handler.boundary_list[:]
                        boundary_list = [1 if q == 0 else q for q in
                                         boundary_list]  # To start at 1
                        color_list = self.handler.color_list[:len(boundary_list) - 1]
                        for b1, b2, c in zip(boundary_list[:-1], boundary_list[1:],
                                             color_list):
                            a.broken_barh([(b1, b2 - b1)],
                                          (low + .25 * (height - low),
                                           .5 * (height - low)),
                                          facecolors=c)
                        fname = 'student {} objective {}.png'.format(user,
                                                                     learn_obj)
                        f.savefig(fname=dirname + fname)
                        matplotlib.pyplot.close()
                        print("coordinates are {}".format(graph_n))
                        print("saved student {} objective {}".format(user,
                                                                     learn_obj))
                        # self.write_spikes(user, learn_obj,
                        #                   self.handler.boundary_list,
                        #                   o_graph, writer)
                    except Exception as e:
                        print("failed saving student {} "
                              "objective {} because of {}".format(user, learn_obj,
                                                                  e))
                        errors += 1
            print("saved all graphs with {} errors".format(errors))

    def write_spikes(self, student, loid, bounds, graph, writer):
        row = [str(student), str(loid)]
        n_peaks, peak_per_bound, trans_peak = self.calc_peaks(graph, bounds)
        print("`1")
        gen_spikiness, spikiness_per_bound = self.calc_spikes(graph, bounds)
        # Spikiness
        row.append(gen_spikiness)
        for sp in spikiness_per_bound:
            row.append(sp)

        # Peaks
        row.append(n_peaks)
        for p in peak_per_bound:
            row.append(p)

        # Transition peaks
        row.append(sum(trans_peak))
        for p in trans_peak:
            row.append(p)

        # Write results
        writer.writerow(row)

    def calc_peaks(self, graph, bounds):
        locs = []
        n_peaks = 0
        peak_per_bound = [0, 0, 0, 0, 0, 0]
        trans_peak = [0, 0, 0, 0, 0, 0]
        m_s = .0015  # Minimum spikiness
        # Track current bound
        bound = 0  # Tracks current bound
        while 0 == bounds[bound+1]:
            bound += 1

        # first answer
        old_j = graph[0]
        if graph[0] > graph[1]+m_s:
            locs.append(0)
            n_peaks += 1
            peak_per_bound[bound] = 1
            trans_peak[bound] = 1

        # middle answers
        for i in range(1, len(graph)-1):
            j = graph[i]
            n = graph[i+1]
            if i+2 > bounds[bound+1]:
                while i+2 > bounds[bound+1]:
                    bound += 1
                    # print("Bound is now at {}".format(bound))
                if (j > old_j+m_s) and (j > n+m_s):
                    # print("Peak, plus transition at position {}".format(i))
                    locs.append(i)
                    trans_peak[bound] += 1
                    n_peaks += 1
                    peak_per_bound[bound] += 1
            else:
                if (j > old_j+m_s) and (j > n+m_s):
                    # print("Peak at position {}".format(i))
                    locs.append(i)
                    n_peaks += 1
                    peak_per_bound[bound] += 1
            old_j = j

        # last answer
        if graph[-1]>old_j+m_s:
            locs.append(len(graph)-1)
            n_peaks+=1
            peak_per_bound[bound] += 1

        print("peaks are at {}".format(locs))
        return n_peaks, peak_per_bound, trans_peak

    def calc_spikes(self, graph, bounds):
        gen_sp = max(graph)/(sum(graph)/len(graph))
        lst_sp = []
        for b1, b2 in zip(bounds[:-1], bounds[1:]):
            part = graph[b1:b2]
            if len(part)==0 or sum(part)==0:
                lst_sp.append("NaN")
            else:
                lst_sp.append(max(part)/(sum(part)/len(part)))
        return gen_sp, lst_sp


class StartPage(tk.Frame):
    """ First page that is shown
	"""

    def __init__(self, parent, controller: GraphGUI):
        tk.Frame.__init__(self, parent)
        self.controller = controller
        # Set up label
        label = ttk.Label(self, text='This is the start page', font=LARGE_FONT)
        label.grid(pady=10, padx=10, columnspan=2, sticky='nswe')

        # Set up menu for user choice.
        options = list(controller.handler.get_users())
        self.variable = tk.StringVar(self)
        self.variable.set(options[0])
        self.user_id = 0

        print(options)
        menu = ttk.OptionMenu(self, self.variable, *options,
                              command=self.controller.upd_user)
        menu.grid(row=1, sticky='we')

        # Set up Button to choose option
        button = ttk.Button(self, text='Get graph for this user',
                            command=lambda: controller.show_frame(GraphPage))
        button.grid(row=2, sticky='we')


class GraphPage(tk.Frame):
    """ First page that is shown

	TODO - Change to real page one
	"""

    def __init__(self, parent, controller):
        tk.Frame.__init__(self, parent)
        # label = ttk.Label(self, text='This is page 1', font=LARGE_FONT)
        # label.pack(pady=10, padx=10)
        # button = ttk.Button(self, text='Go back home',
        # 	command=lambda: controller.show_frame(StartPage))
        # button.pack()
        self.controller = controller

        # Set up menu for user choice.
        stud_id_options = list(controller.handler.get_users())
        self.stud_id_variable = tk.StringVar(self)
        self.stud_id_variable.set(controller.user_id)

        learn_obj_id_options = \
            list(np.unique(controller.handler.learn_obj_ids)[1:])
        self.learn_obj_id_variable = tk.StringVar(self)
        # self.learn_obj_id_variable.set(learn_obj_id_options[0])

        navigation_frame = tk.Frame(self)
        button1 = ttk.Button(navigation_frame, text='<',
                             command=lambda: self.next_user(True))
        menu = ttk.OptionMenu(navigation_frame, self.stud_id_variable,
                              stud_id_options[0], *stud_id_options,
                              command=controller.upd_user)
        learning_goal_menu = ttk.OptionMenu(navigation_frame,
                                            self.learn_obj_id_variable,
                                            learn_obj_id_options[0],
                                            *learn_obj_id_options,
                                            command=controller.upd_learn_goal)
        button2 = ttk.Button(navigation_frame, text='>',
                             command=lambda: self.next_user())
        button1.grid(column=0, row=0)
        ttk.Label(navigation_frame, text="Pick a student:").grid(column=1,
                                                                 row=0)
        menu.grid(column=2, row=0)
        button2.grid(column=3, row=0)
        ttk.Label(navigation_frame, text="Pick the learning goal:").grid(
            column=4,
            row=0)
        learning_goal_menu.grid(column=5, row=0)
        navigation_frame.pack()

        method_frame = tk.Frame(self)
        ttk.Label(method_frame, text='Include attempts of exercise:').grid(
            row=0,
            column=0)
        self.method_var = tk.IntVar()
        self.method_var.set(0)
        self.methods = ['all', 'first', 'second', 'all but first', 'last']
        for i, pick in enumerate(self.methods):
            chk = ttk.Radiobutton(method_frame, text=pick,
                                  variable=self.method_var,
                                  value=i, command=self.change_method)
            chk.grid(row=0, column=i + 1)
        ttk.Button(method_frame, text='Save graph',
                   command=controller.save_graph).grid(row=0, column=len(
            self.methods) + 1)
        ttk.Button(method_frame,
                   text='Save all graphs',
                   command=controller.save_all_graphs).grid(row=0,
                                                            column=len(
                                                                self.methods) + 2)
        method_frame.pack()

        canvas = FigureCanvasTkAgg(controller.f, self)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)
        canvas._tkcanvas.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        values_frame = tk.Frame(self)
        label_l0 = ttk.Label(values_frame, text='L0:')
        label_t = ttk.Label(values_frame, text='T:')
        label_g = ttk.Label(values_frame, text='G:')
        label_s = ttk.Label(values_frame, text='S:')
        self.l0_variable = tk.StringVar()
        self.l0_variable.set(str(controller.handler.m2m.p_l0))
        textbox_l0 = ttk.Entry(values_frame, textvariable=self.l0_variable)
        self.t_variable = tk.StringVar()
        self.t_variable.set(str(controller.handler.m2m.p_T))
        textbox_t = ttk.Entry(values_frame, textvariable=self.t_variable)
        self.g_variable = tk.StringVar()
        self.g_variable.set(str(controller.handler.m2m.p_G))
        textbox_g = ttk.Entry(values_frame, textvariable=self.g_variable)
        self.s_variable = tk.StringVar()
        self.s_variable.set(str(controller.handler.m2m.p_S))
        textbox_s = ttk.Entry(values_frame, textvariable=self.s_variable)
        label_l0.grid(row=0, column=0)
        textbox_l0.grid(row=0, column=1)
        label_t.grid(row=0, column=2)
        textbox_t.grid(row=0, column=3)
        label_g.grid(row=0, column=4)
        textbox_g.grid(row=0, column=5)
        label_s.grid(row=0, column=6)
        textbox_s.grid(row=0, column=7)
        button_update = ttk.Button(values_frame, text='Update values',
                                   command=self.update_values)
        button_update.grid(row=0, column=8)
        values_frame.pack()

    def update_values(self):
        m2m_list = [self.controller.handler.m2m.p_l0,
                    self.controller.handler.m2m.p_T,
                    self.controller.handler.m2m.p_G,
                    self.controller.handler.m2m.p_S]
        for id, entry in enumerate([self.l0_variable, self.t_variable,
                                    self.g_variable, self.s_variable]):
            try:
                value = float(entry.get())
            except ValueError:
                value = 0.00
                entry.set(0.00)
            m2m_list[id] = value
            entry.set(value)
            print(m2m_list)
        self.controller.handler.m2m.p_l0 = m2m_list[0]
        self.controller.handler.m2m.p_T = m2m_list[1]
        self.controller.handler.m2m.p_G = m2m_list[2]
        self.controller.handler.m2m.p_S = m2m_list[3]
        self.controller.update_graph(0)

    def change_method(self):
        method = self.methods[self.method_var.get()]
        self.controller.method = method

    def next_user(self, reverse=False):
        id_plus = 1
        if reverse is True:
            id_plus = -1
        self.stud_id_variable.set(self.controller.move_user(id_plus))


def qf(quick_print):
    print(quick_print)


if __name__ == '__main__':
    app = GraphGUI()
    ani = animation.FuncAnimation(app.f, app.update_graph, interval=100)
    app.mainloop()
