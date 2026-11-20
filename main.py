"""Interactive double pendulum. Run: python main.py"""
import sys
from PyQt6.QtWidgets import QApplication
from visualization import MainWindow


def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == '__main__':
    sys.exit(main())
