def length_to_string(length): 
    if length < 1000: 
        return f'{length}bp'
    elif length < 1000000: 
        return f'{int(length/1000)}kb'
    else: 
        return f'{int(length/1000000)}Mb'
    
