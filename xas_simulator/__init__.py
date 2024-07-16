from . import xas 

def plugin_menu():
    menu = 'XAS Simulator'
    actions = []
    actions.append(('XAS setup', xas.show_dialog))
    #actions.append(('Convert to Q-E', convert_qe.show_dialog))
    return menu, actions