import turtle
import os
from PIL import Image 

scale = 4

Lport = 105*scale
Lhalf = 50*scale
Lslit = Lport - 2*Lhalf
Lcase = 118.9*scale
Dport = 25*scale
Dout = 76.8*scale
Dt = 8.37*scale
din = 1.5*scale

a = 600
b = 400
size = 3
speed = 0

bgcolor = "white"
pencolor_case = "black"
fillcolor_case = "silver"
pencolor_grain = "orange"
fillcolor_grain = "moccasin"
pencolor_inhibitor = "crimson"
fillcolor_inhibitor = "coral"

def draw1(R1, R2, L, stay = True):

    window = turtle.Screen()
    window.screensize(a, b)
    window.bgcolor(bgcolor)

    pen = turtle.Turtle()
    pen.hideturtle()
    pen.pensize(size)
    pen.speed(speed)

    # Draw Case
    draw_case(pen)

    # Draw Grain
    draw_grain(pen, R1, R2, L)

    # Draw Inhibitor
    draw_inhibitor(pen)

    if stay:
        turtle.done()

    return window

def draw2(R11, R12, R21, R22, L1, L2, stay = True):
    
    window = turtle.Screen()
    window.screensize(a, b)
    window.bgcolor(bgcolor)

    pen = turtle.Turtle()
    pen.hideturtle()
    pen.pensize(size)
    pen.speed(speed)

    # Draw Case
    draw_case(pen)

    # Draw Grain
    draw_grain1(pen, R11, R12, L1)
    draw_grain2(pen, R21, R22, L2)

    if stay:
        turtle.done()

    return window

def draw_case(pen):
    pen.color(pencolor_case)
    pen.fillcolor(fillcolor_case)

    pen.penup()
    pen.goto(Lcase/2 + size, Dt/2 + size)
    pen.pendown()

    pen.begin_fill()
    pen.goto(Lcase/2 + size, Dout/2 + size)
    pen.goto(-Lcase/2 - size, Dout/2 + size)
    pen.goto(-Lcase/2 - size, -Dout/2 - size)
    pen.goto(Lcase/2 + size, -Dout/2 - size)
    pen.goto(Lcase/2 + size, -Dt/2 - size)
    pen.end_fill()

    turtle.update()

def draw_grain(pen, R1, R2, L):
    pen.color(pencolor_grain)
    pen.fillcolor(fillcolor_grain)

    # Top
    pen.penup()
    pen.goto(-Lport/2, R1)
    pen.pendown()
    
    pen.begin_fill()
    pen.goto(-Lport/2, Dout/2)
    pen.goto(-Lport/2 + L, Dout/2)
    pen.goto(-Lport/2 + L, R2)
    pen.goto(-Lport/2, R1)
    pen.end_fill()

    # Bottom
    pen.penup()
    pen.goto(-Lport/2, -R1)
    pen.pendown()
    
    pen.begin_fill()
    pen.goto(-Lport/2, -Dout/2)
    pen.goto(-Lport/2 + L, -Dout/2)
    pen.goto(-Lport/2 + L, -R2)
    pen.goto(-Lport/2, -R1)
    pen.end_fill()

    turtle.update()

def draw_grain1(pen, R1, R2, L):
    pen.color(pencolor_grain)
    pen.fillcolor(fillcolor_grain)

    # Top
    pen.penup()
    pen.goto(-Lslit/2 - Lhalf/2 - L/2, R1)
    pen.pendown()
    
    pen.begin_fill()
    pen.goto(-Lslit/2 - Lhalf/2 - L/2, Dout/2)
    pen.goto(-Lslit/2 - Lhalf/2 + L/2, Dout/2)
    pen.goto(-Lslit/2 - Lhalf/2 + L/2, R2)
    pen.goto(-Lslit/2 - Lhalf/2 - L/2, R1)
    pen.end_fill()

    # Bottom
    pen.penup()
    pen.goto(-Lslit/2 - Lhalf/2 - L/2, -R1)
    pen.pendown()
    
    pen.begin_fill()
    pen.goto(-Lslit/2 - Lhalf/2 - L/2, -Dout/2)
    pen.goto(-Lslit/2 - Lhalf/2 + L/2, -Dout/2)
    pen.goto(-Lslit/2 - Lhalf/2 + L/2, -R2)
    pen.goto(-Lslit/2 - Lhalf/2 - L/2, -R1)
    pen.end_fill()

    turtle.update()

def draw_grain2(pen, R1, R2, L):
    pen.color(pencolor_grain)
    pen.fillcolor(fillcolor_grain)

    # Top
    pen.penup()
    pen.goto(Lslit/2 + Lhalf/2 - L/2, R1)
    pen.pendown()
    
    pen.begin_fill()
    pen.goto(Lslit/2 + Lhalf/2 - L/2, Dout/2)
    pen.goto(Lslit/2 + Lhalf/2 + L/2, Dout/2)
    pen.goto(Lslit/2 + Lhalf/2 + L/2, R2)
    pen.goto(Lslit/2 + Lhalf/2 - L/2, R1)
    pen.end_fill()

    # Bottom
    pen.penup()
    pen.goto(Lslit/2 + Lhalf/2 - L/2, -R1)
    pen.pendown()
    
    pen.begin_fill()
    pen.goto(Lslit/2 + Lhalf/2 - L/2, -Dout/2)
    pen.goto(Lslit/2 + Lhalf/2 + L/2, -Dout/2)
    pen.goto(Lslit/2 + Lhalf/2 + L/2, -R2)
    pen.goto(Lslit/2 + Lhalf/2 - L/2, -R1)
    pen.end_fill()

    turtle.update()

def draw_inhibitor(pen):
    pen.color(pencolor_inhibitor)
    pen.fillcolor(fillcolor_inhibitor)

    # Top Left
    pen.penup()
    pen.goto(-Lport/2 - size, Dport/2)
    pen.pendown()

    pen.begin_fill()
    pen.goto(-Lport/2 - size - din - size, Dport/2)
    pen.goto(-Lport/2 - size - din - size, Dout/2)
    pen.goto(-Lport/2 - size, Dout/2)
    pen.goto(-Lport/2 - size, Dport/2)
    pen.end_fill()

    # Top Right
    pen.penup()
    pen.goto(Lport/2 + size, Dport/2)
    pen.pendown()

    pen.begin_fill()
    pen.goto(Lport/2 + size + din + size, Dport/2)
    pen.goto(Lport/2 + size + din + size, Dout/2)
    pen.goto(Lport/2 + size, Dout/2)
    pen.goto(Lport/2 + size, Dport/2)
    pen.end_fill()

    # Bottom Left
    pen.penup()
    pen.goto(-Lport/2 - size, -Dport/2)
    pen.pendown()

    pen.begin_fill()
    pen.goto(-Lport/2 - size - din - size, -Dport/2)
    pen.goto(-Lport/2 - size - din - size, -Dout/2)
    pen.goto(-Lport/2 - size, -Dout/2)
    pen.goto(-Lport/2 - size, -Dport/2)
    pen.end_fill()

    # Bottom Right
    pen.penup()
    pen.goto(Lport/2 + size, -Dport/2)
    pen.pendown()

    pen.begin_fill()
    pen.goto(Lport/2 + size + din + size, -Dport/2)
    pen.goto(Lport/2 + size + din + size, -Dout/2)
    pen.goto(Lport/2 + size, -Dout/2)
    pen.goto(Lport/2 + size, -Dport/2)
    pen.end_fill()

    turtle.update()

def frames1(R1, R2, L):
    for i in range(len(R1)):
        window = draw1(R1[i], R2[i], L[i], stay = False)
        fileName = "Frame01%03d" % (i)

        # Save as .eps
        window.setup(a, b)
        window.tracer(False)
        window.tracer(True)
        canvas = window.getcanvas()
        canvas.postscript(file = "Frames/eps/" + fileName + ".eps", width = a, height = b)
        
        # Convert to .jpg
        img = Image.open("Frames/eps/" + fileName + ".eps") 
        img.load(scale = 10)
        img.save("Frames/jpg/" + fileName + ".jpg")

def frames2(R11, R12, R21, R22, L1, L2):
    for i in range(len(R11)):
        window = draw2(R11[i], R12[i], R21[i], R22[i], L1[i], L2[i], stay = False)
        fileName = "Frame02%03d" % (i)

        # Save as .eps
        window.setup(a, b)
        window.tracer(False)
        window.tracer(True)
        canvas = window.getcanvas()
        canvas.postscript(file = "Frames/eps/" + fileName + ".eps", width = a, height = b)
        
        # Convert to .jpg
        img = Image.open("Frames/eps/" + fileName + ".eps") 
        img.load(scale = 10)
        img.save("Frames/jpg/" + fileName + ".jpg")

R1 = [1*Dport/2, 1.1*Dport/2, 1.2*Dport/2, 1.3*Dport/2]
R2 = [1*Dport/2, 1.2*Dport/2, 1.4*Dport/2, 1.6*Dport/2]
#draw(R1[2], R2[2])
#frames1(R1, R2)

def animate1(R1, R2, L):
    R1 *= 1000*scale
    R2 *= 1000*scale
    L *= 1000*scale

    frames1(R1, R2, L)

    os.system("ffmpeg -f image2 -r 1/0.1 -i ./Frames/jpg/Frame01%03d.jpg -vcodec mpeg4 -y ./Configuration1.mp4")

    return 1

def animate2(R11, R12, R21, R22, L1, L2):
    R11 *= 1000*scale
    R12 *= 1000*scale
    R21 *= 1000*scale
    R22 *= 1000*scale
    L1 *= 1000*scale
    L2 *= 1000*scale

    frames2(R11, R12, R21, R22, L1, L2)

    os.system("ffmpeg -f image2 -r 1/0.1 -i ./Frames/jpg/Frame02%03d.jpg -vcodec mpeg4 -y ./Configuration2.mp4")

    return 1
