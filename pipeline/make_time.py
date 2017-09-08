from time import sleep
i = 0
while i < 3600:
    i = i + 1
    sleep(1)
    print('remaining time {}'.format(3600-i))

