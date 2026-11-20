"""PyQt6 interface, interactive pendulum and live pyqtgraph plots."""
import csv

import numpy as np
import pyqtgraph as pg
from PyQt6.QtCore import Qt, QTimer, QPointF, pyqtSignal
from PyQt6.QtGui import QPainter, QPen, QBrush, QColor
from PyQt6.QtWidgets import (QMainWindow, QWidget, QHBoxLayout, QVBoxLayout,
    QFormLayout, QDoubleSpinBox, QPushButton, QCheckBox, QLabel,
    QFileDialog, QMessageBox, QGroupBox, QScrollArea, QGridLayout)

from controller import PIDController, suggested_gains, upper_linearization
from integrator import rk4_step
from physics import derivatives, energy_scale, total_energy

DEFAULTS = dict(m1=1., m2=1., L1=1., L2=1., g=9.81, c1=0., c2=0.,
                theta1_deg=120., theta2_deg=-35., omega1_deg=0., omega2_deg=0., dt=0.005)


class PendulumView(QWidget):
    dragged = pyqtSignal(int, float)

    def __init__(self):
        super().__init__()
        self.setMinimumSize(330, 350)
        self.state = np.radians([120., 0., -35., 0.])
        self.lengths = (1., 1.)
        self.drag_index = None
        self.setToolTip('В остановленном состоянии и на паузе перетаскивайте грузы мышью')

    def geometry_points(self):
        w, h = self.width(), self.height()
        scale = 0.43 * min(w, h) / max(sum(self.lengths), 0.01)
        center = QPointF(w / 2, h / 2)
        a, b = self.state[0], self.state[2]
        p1 = QPointF(center.x() + scale*self.lengths[0]*np.sin(a), center.y() + scale*self.lengths[0]*np.cos(a))
        p2 = QPointF(p1.x() + scale*self.lengths[1]*np.sin(b), p1.y() + scale*self.lengths[1]*np.cos(b))
        return center, p1, p2

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.fillRect(self.rect(), QColor('#fafcff'))
        p0, p1, p2 = self.geometry_points()
        p.setPen(QPen(QColor('#b4c1cf'), 1, Qt.PenStyle.DashLine))
        p.drawLine(int(p0.x()), 0, int(p0.x()), self.height())
        p.setPen(QPen(QColor('#334155'), 4))
        p.drawLine(p0, p1)
        p.drawLine(p1, p2)
        for point, color, radius in ((p0, '#334155', 5), (p1, '#0284c7', 9), (p2, '#ef4444', 10)):
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QBrush(QColor(color)))
            p.drawEllipse(point, radius, radius)
        p.end()

    def mousePressEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        _, p1, p2 = self.geometry_points()
        pos = event.position()
        if (pos-p2).manhattanLength() < 32:
            self.drag_index = 2
        elif (pos-p1).manhattanLength() < 32:
            self.drag_index = 1

    def mouseMoveEvent(self, event):
        if self.drag_index is None:
            return
        p0, p1, _ = self.geometry_points()
        anchor = p0 if self.drag_index == 1 else p1
        delta = event.position() - anchor
        angle = float(np.degrees(np.arctan2(delta.x(), delta.y())))
        self.dragged.emit(self.drag_index, angle)

    def mouseReleaseEvent(self, event):
        self.drag_index = None


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Двойной маятник — PyQt6 / RK4 / PID')
        self.resize(1450, 900)
        self.timer = QTimer(self)
        self.timer.setInterval(20)
        self.timer.timeout.connect(self.tick)
        self.running = False
        self.elapsed = 0.
        self.records = []
        self.max_plot_points = 2500
        self.view = PendulumView()
        self.controls = {}
        self.build_ui()
        self.pid_enabled.setChecked(False)
        self.autotune()
        self.reset()

    def spin(self, name, value, low, high, decimals=3, step=0.1):
        s = QDoubleSpinBox()
        s.setRange(low, high)
        s.setDecimals(decimals)
        s.setSingleStep(step)
        s.setValue(value)
        self.controls[name] = s
        return s

    def build_ui(self):
        root = QWidget()
        self.setCentralWidget(root)
        layout = QHBoxLayout(root)
        panel = QWidget()
        form_layout = QVBoxLayout(panel)
        form_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        phys = QGroupBox('Физика и начальные условия')
        f = QFormLayout(phys)
        labels = {'m1':'Масса 1, кг', 'm2':'Масса 2, кг', 'L1':'Длина 1, м',
                  'L2':'Длина 2, м', 'g':'Гравитация, м/с²', 'c1':'Трение опоры',
                  'c2':'Трение шарнира', 'theta1_deg':'Угол 1, °', 'theta2_deg':'Угол 2, °',
                  'omega1_deg':'Скорость 1, °/с', 'omega2_deg':'Скорость 2, °/с', 'dt':'Шаг RK4, с'}
        for key, val in DEFAULTS.items():
            low, high = ((0.01, 1000.) if key in ('m1','m2','L1','L2','g') else
                         (0., 100.) if key in ('c1','c2') else
                         (0.0001, 0.05) if key == 'dt' else (-3600., 3600.))
            f.addRow(labels[key], self.spin(key, val, low, high, 4 if key == 'dt' else 2, 0.001 if key == 'dt' else 0.1))
        form_layout.addWidget(phys)
        control_box = QGroupBox('PID: стабилизация верхнего положения')
        pf = QFormLayout(control_box)
        self.pid_enabled = QCheckBox('Включить PID (ТОЛЬКО возле 180°, 180°)')
        pf.addRow(self.pid_enabled)
        self.pid_warning = QLabel('По умолчанию свободное движение. PID действует только при отклонении каждого звена не более 25° и скорости не более 90°/с; при выходе за пределы отключается. Для стабилизации задайте углы около 180°.')
        self.pid_warning.setWordWrap(True)
        pf.addRow(self.pid_warning)
        self.gain_spins = []
        for name in ('P₁ (угол 1)', 'D₁ (скорость 1)', 'P₂ (угол 2)', 'D₂ (скорость 2)'):
            s = QDoubleSpinBox()
            s.setRange(-100000, 100000)
            s.setDecimals(3)
            self.gain_spins.append(s)
            pf.addRow(name, s)
        self.ki_spin = self.spin('ki', 0., -10000., 10000.)
        self.limit_spin = self.spin('limit', 100., 0.01, 100000.)
        pf.addRow('Ki (интеграл угла 1)', self.ki_spin)
        pf.addRow('Макс. момент, Н·м', self.limit_spin)
        auto = QPushButton('Подобрать PD по полюсам')
        auto.clicked.connect(self.autotune)
        pf.addRow(auto)
        form_layout.addWidget(control_box)
        buttons = QGridLayout()
        for i, (label, handler) in enumerate((('Старт', self.start), ('Пауза', self.pause),
            ('Стоп / сброс', self.reset), ('Save to PNG', self.save_png), ('Save to CSV', self.save_csv))):
            btn = QPushButton(label)
            btn.clicked.connect(handler)
            buttons.addWidget(btn, i//2, i%2)
        form_layout.addLayout(buttons)
        self.status = QLabel('Готово')
        self.status.setWordWrap(True)
        form_layout.addWidget(self.status)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(panel)
        scroll.setMinimumWidth(320)
        scroll.setMaximumWidth(385)
        layout.addWidget(scroll)
        right = QWidget()
        rlayout = QVBoxLayout(right)
        rlayout.addWidget(self.view, 4)
        self.view.dragged.connect(self.on_drag)
        plots = QGridLayout()
        pg.setConfigOptions(antialias=True, background='w', foreground='k')
        self.plots = []
        titles = ('Углы от времени', 'Изменение энергии (при трении/PID не сохраняется)', 'Фазовый портрет 1', 'Фазовый портрет 2')
        for i, title in enumerate(titles):
            plot = pg.PlotWidget(title=title)
            plot.showGrid(x=True, y=True, alpha=0.25)
            self.plots.append(plot)
            plots.addWidget(plot, i//2, i%2)
        self.angle_lines = [self.plots[0].plot(pen=pg.mkPen(c, width=2)) for c in ('#0284c7', '#ef4444')]
        self.energy_line = self.plots[1].plot(pen=pg.mkPen('#334155', width=2))
        self.phase_lines = [self.plots[i].plot(pen=pg.mkPen('#0284c7', width=1.5)) for i in (2,3)]
        self.plots[0].setLabel('bottom', 't', 'с')
        self.plots[0].setLabel('left', 'угол', '°')
        self.plots[1].setLabel('bottom', 't', 'с')
        for i in (2,3):
            self.plots[i].setLabel('bottom', 'угол', '°')
            self.plots[i].setLabel('left', 'скорость', '°/с')
        rlayout.addLayout(plots, 5)
        layout.addWidget(right, 1)
        for name in ('theta1_deg', 'theta2_deg'):
            self.controls[name].valueChanged.connect(self.preview)
        for name in ('L1', 'L2'):
            self.controls[name].valueChanged.connect(self.preview)

    def parameters(self):
        return {k: v.value() for k,v in self.controls.items() if k in DEFAULTS}

    def preview(self):
        if self.running or self.elapsed > 0:
            return
        p = self.parameters()
        self.view.lengths = (p['L1'], p['L2'])
        self.view.state = np.radians([p['theta1_deg'], p['omega1_deg'], p['theta2_deg'], p['omega2_deg']])
        self.view.update()

    def on_drag(self, index, angle):
        if self.running:
            return
        if self.elapsed > 0.:
            # On pause, explicitly set current geometry and stop angular speed.
            j = 0 if index == 1 else 2
            self.state[j] = np.deg2rad(angle)
            self.state[j + 1] = 0.0
            if self.controller:
                self.controller.reset()
            self.view.state = self.state.copy()
            self.view.update()
            self.append_record(0.)
            self.redraw()
        else:
            self.controls['theta1_deg' if index == 1 else 'theta2_deg'].setValue(angle)


    def autotune(self):
        try:
            gains = suggested_gains(self.parameters())
            for spin, value in zip(self.gain_spins, gains):
                spin.setValue(float(value))
            self.status.setText('PD коэффициенты рассчитаны размещением полюсов верхней линеаризации (не LQR).')
        except Exception as exc:
            QMessageBox.warning(self, 'Настройка PID', str(exc))

    def reset(self):
        self.timer.stop()
        self.running = False
        self.elapsed = 0.
        self.p = self.parameters()
        self.state = np.radians([self.p['theta1_deg'], self.p['omega1_deg'], self.p['theta2_deg'], self.p['omega2_deg']])
        self.controller = None
        self.records = []
        self.append_record(0.)
        self.preview()
        self.redraw()
        self.status.setText('Остановлено / сброшено. Можно менять параметры и перетаскивать грузы.')

    def start(self):
        if self.running:
            return
        if self.elapsed == 0.:
            self.p = self.parameters()
            self.state = np.radians([self.p['theta1_deg'], self.p['omega1_deg'], self.p['theta2_deg'], self.p['omega2_deg']])
            self.records = []
            self.controller = PIDController([s.value() for s in self.gain_spins], self.ki_spin.value(), self.limit_spin.value()) if self.pid_enabled.isChecked() else None
            if self.controller and not self.controller.is_local(self.state):
                self.status.setText('Начальная поза вне области PID: управление отключено, будет свободное движение.')
            self.append_record(0.)
        self.running = True
        self.timer.start()
        self.status.setText('Запущено: свободное движение' if self.controller is None else 'Запущено: локальный PID (вне области автоматически выключается)')

    def pause(self):
        self.timer.stop()
        self.running = False
        self.status.setText('Пауза. Старт продолжит симуляцию; Стоп сбросит.')

    def append_record(self, tau):
        E = float(total_energy(self.state, self.p)[0])
        self.records.append((self.elapsed, *self.state, tau, E))

    def tick(self):
        try:
            dt = self.p['dt']
            # Fixed number of RK4 steps per GUI frame; physics is independent of frame rate.
            n = max(1, int(round(0.02/dt)))
            for _ in range(n):
                tau = self.controller.output(self.state) if self.controller else 0.
                self.state = rk4_step(lambda y, t: derivatives(y, t, self.p, tau=tau), self.state, self.elapsed, dt)
                if not np.all(np.isfinite(self.state)):
                    raise FloatingPointError('Нефинитное состояние: уменьшите шаг интегрирования.')
                if self.controller:
                    self.controller.advance(self.state, dt)
                self.elapsed += dt
                self.append_record(tau)
            self.view.state = self.state.copy()
            self.view.lengths = (self.p['L1'], self.p['L2'])
            self.view.update()
            self.redraw()
            mode = ('Свободное движение' if self.controller is None else
                    'PID отключён (вне области)' if self.controller.disengaged else 'PID')
            self.status.setText(f'{mode} | t = {self.elapsed:.2f} с | момент = {tau:.3f} Н·м | точек CSV: {len(self.records)}')
        except Exception as exc:
            self.pause()
            QMessageBox.critical(self, 'Ошибка симуляции', str(exc))

    def redraw(self):
        # Keep full history for CSV; downsample only the plots for UI responsiveness.
        data = np.asarray(self.records, dtype=float)
        step = max(1, len(data)//self.max_plot_points)
        d = data[::step]
        t = d[:,0]
        angles = np.rad2deg(d[:,[1,3]])
        for i in range(2):
            self.angle_lines[i].setData(t, angles[:,i])
        scale = energy_scale(self.p)
        self.energy_line.setData(t, (d[:,6]-data[0,6])/scale)
        for i, (a,w) in enumerate(((1,2),(3,4))):
            self.phase_lines[i].setData(np.rad2deg(d[:,a]), np.rad2deg(d[:,w]))

    def save_png(self):
        path, _ = QFileDialog.getSaveFileName(self, 'Сохранить PNG', 'pendulum.png', 'PNG (*.png)')
        if not path:
            return
        from PyQt6.QtGui import QPixmap
        target = self.centralWidget()
        pix = target.grab()
        if not pix.save(path, 'PNG'):
            QMessageBox.warning(self, 'Ошибка', 'Не удалось сохранить PNG')

    def save_csv(self):
        path, _ = QFileDialog.getSaveFileName(self, 'Сохранить CSV', 'pendulum.csv', 'CSV (*.csv)')
        if not path:
            return
        with open(path, 'w', newline='', encoding='utf-8-sig') as f:
            writer = csv.writer(f)
            writer.writerow(['time_s','theta1_rad','omega1_rad_s','theta2_rad','omega2_rad_s','torque_Nm','energy_J'])
            writer.writerows(self.records)
