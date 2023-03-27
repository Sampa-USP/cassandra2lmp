def write_last(l_frame,ofile = "/tmp/h2o.pdb"):
    lista = ['','','','']
    count = 0
    for indice,atom in enumerate(l_frame):
        _symbol = atom.symbol
        lista[count] = _symbol
        
        #print(_symbol)
        count +=1
        count = count % 4
        
        if lista == ['O','H','H','O']:
            break

    sio2 = l_frame[:indice - 3]
    h2o = l_frame[indice - 3:]

    h2o.write(ofile)