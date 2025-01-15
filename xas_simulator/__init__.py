from . import xas
from . import xmcd_loader
from . import nexpy_nb_runner
# from . import peem_loader


def plugin_menu():
    menu = 'XAS Simulator'
    actions = []
    actions.append(('XAS setup', xas.show_dialog))
    actions.append(('Data loader', xmcd_loader.show_dialog))
    actions.append(('Notebook Runner', nexpy_nb_runner.show_dialog))
    # actions.append(('PEEM loader', peem_loader.show_dialog))
    #actions.append(('Convert to Q-E', convert_qe.show_dialog))
    return menu, actions
