from datetime import datetime

def vprint(*args, **kwargs):
    vstr = (datetime.now().strftime("%m-%d-%Y %H:%M:%S.%f"))[:-3]
    print("[VorAMR {}] ".format(vstr), end='')
    return print(*args, **kwargs)
