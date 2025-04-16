def fun(b,*args,a ):
    print("a={}".format(a))
    print("b={}".format(b))
    for arg in args:
        print('Optional argument: {}'.format( arg ) )

fun(24,[1,22],[2,3],a=33)