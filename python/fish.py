__FISH__ = '''o
o      ______/~/~/~/__           /((
  o  // __            ====__    /_((
 o  //  @))       ))))      ===/__((
    ))           )))))))        __((
    \\     \)     ))))    __===\ _((
     \\_______________====      \_((
                                 \(('''
                                 
def blub(msg):
    h = "_" * len(msg)
    e = " " * len(msg)
    print(f" _{h}_ ")
    #print(f"/ {e} \\")
    print(f"/ {msg} \\")
    print(f"\\_{h}_/")
    print(__FISH__)

if __name__ == '__main__':
  blub("long test message")