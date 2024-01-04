class Messages:
    @staticmethod
    def firstCorona():
        message = 'We are now computing the 1st corona'
        print(message)
    
    @staticmethod
    def laterCorona(index):
        message = f'We are now computing the {index+2}nd corona'
        print(message)

    @staticmethod
    def heeschNumber(index):
        message = f'The Heesch number is {index + 1}'
        print(message)

    @staticmethod
    def coronaCount(currentIndex, totalAmount):
        message = f'We are at corona {currentIndex} out of {totalAmount}\r'
        print( message, end='' )

    @staticmethod
    def totalNewCoronaAmount(totalNewAmount):
        messages = f'In total there are {totalNewAmount} working configurations'
        print(messages)