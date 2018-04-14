def info(object, docString, omitFunctionName=0, spacing=0, collapse=0, blanklines=0): 
    """Print methods and doc strings.
    Takes module, class, list, dictionary, or string."""
    
    methodList = [method for method in dir(object) if callable(getattr(object, method))]
    
    nameFunc = omitFunctionName and (lambda s: "") or (lambda s: s) 
    if collapse:
        docFunc = (lambda s: " ".join(s.split()))
    else:
        docFunc = (lambda s: s.replace("\n", "\n   ")) 
    lineFunc = (blanklines and (lambda: "\n")) or (lambda: "")

    print "\n%s:%s" % (docString, lineFunc())
    print "\n".join([" * %s%s%s" %
                    (nameFunc(method.ljust(spacing)+"() "),
                     docFunc(str(getattr(object, method).__doc__)),
                     lineFunc())
                     for method in methodList]) 
