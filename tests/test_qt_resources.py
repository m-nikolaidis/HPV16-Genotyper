import os
import unittest


os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtGui import QPixmap
from PyQt5.QtCore import QSize
from PyQt5.QtWidgets import QApplication, QMainWindow, QVBoxLayout, QPushButton, QWidget

from hpv16genotyper import app


class QtResourceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.qt_app = QApplication.instance() or QApplication([])

    def test_magnifying_glass_icon_is_registered(self):
        self.assertFalse(QPixmap(":/resources/icons/cil-magnifying-glass.png").isNull())

    def test_menu_icon_uses_fixed_size_qicon(self):
        button = QPushButton()
        app.GuiFunctions._setMenuIcon(
            button, "url(:/resources/icons/cil-home.png)"
        )

        self.assertFalse(button.icon().isNull())
        self.assertEqual(button.iconSize(), QSize(24, 24))
        self.assertNotIn("background-image", app.Style.style_bt_standard)

    def test_main_page_menu_receives_home_icon(self):
        class MenuHost(QWidget):
            def __init__(self):
                super().__init__()
                self.buttons = {}
                self.menusLayout = QVBoxLayout(self)
                self.ui = type("Ui", (), {"menusLayout": self.menusLayout})()

            def Button(self):
                pass

        host = MenuHost()
        app.GuiFunctions.addNewMenu(
            host,
            "Main Page",
            "homeButton",
            ":/resources/icons/cil-home.png",
            True,
        )

        button = host.buttons["Main Page"]
        self.assertFalse(button.icon().isNull())
        self.assertEqual(button.iconSize(), QSize(24, 24))
        self.assertEqual(button.text(), "")
        self.assertEqual(button.minimumWidth(), 70)
        self.assertNotIn("padding-left: 45px", app.Style.style_sidebar_menu)

    def test_left_menu_keeps_compact_width(self):
        window = QMainWindow()
        ui = app.Ui_MainWindow()
        ui.setupUi(window)

        self.assertEqual(ui.frame_left_menu.minimumWidth(), 70)
        self.assertEqual(ui.frame_left_menu.maximumWidth(), 70)


if __name__ == "__main__":
    unittest.main()
