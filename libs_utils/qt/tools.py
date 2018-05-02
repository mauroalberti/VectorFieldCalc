
from PyQt5.QtWidgets import *


def info(parent, header, msg):
    """

    :param parent:
    :param header:
    :param msg:
    :return:
    """
    
    QMessageBox.information(parent, header, msg)


def warn(parent, header, msg):
    """

    :param parent:
    :param header:
    :param msg:
    :return:
    """

    QMessageBox.warning(parent, header, msg)


def error(parent, header, msg):
    """

    :param parent:
    :param header:
    :param msg:
    :return:
    """

    QMessageBox.error(parent, header, msg)
    
    
def update_ComboBox(combobox, init_choice, names):
    """

    :param combobox:
    :param init_choice:
    :param names:
    :return:
    """

    combobox.clear()

    if len(names) == 0:
        return

    if init_choice:
        combobox.addItem(init_choice)

    combobox.addItems(names)

