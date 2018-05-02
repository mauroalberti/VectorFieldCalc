# -*- coding: utf-8 -*-

from PyQt5.QtCore import *
from PyQt5.QtWidgets import *


def lastUsedDir(module_nm):
    """
    from module RASTERCALC by Barry Rowlingson

    :return:
    """

    settings = QSettings()
    return settings.value("/{}/lastDir".format(module_nm), "", type=str)


def setLastUsedDir(module_nm, lastDir):
    """
    from module RASTERCALC by Barry Rowlingson

    :return:
    """

    path = QFileInfo(lastDir).absolutePath()
    settings = QSettings()
    settings.setValue("/{}/lastDir".format(module_nm), str(path))


def update_directory_key(settings, settings_dir_key, fileName):
    """
    modified from module RASTERCALC by Barry Rowlingson
    """

    path = QFileInfo(fileName).absolutePath()
    settings.setValue(settings_dir_key,
                      str(path))


def new_file_path(parent, show_msg, path, filter_text):
    """

    :param parent:
    :param show_msg:
    :param path:
    :param filter_text:
    :return:
    """

    output_filename = QFileDialog.getSaveFileName(parent,
                                                  show_msg,
                                                  path,
                                                  filter_text)
    if not output_filename:
        return ''
    else:
        return output_filename


def old_file_path(parent, show_msg, filter_extension, filter_text):
    """

    :param parent:
    :param show_msg:
    :param filter_extension:
    :param filter_text:
    :return:
    """

    input_filename = QFileDialog.getOpenFileName(parent,
                                                 parent.tr(show_msg),
                                                 filter_extension,
                                                 filter_text)
    if not input_filename:
        return ''
    else:
        return input_filename

