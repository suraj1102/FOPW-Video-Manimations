"""
Manim Animations explaining Copenhagen Interpretation and 
Many Worlds Theory 

Created by Suraj Dayma - github.com/suraj1102
suraj.dayma.ug23@plaksha.edu.in

References Used:
Manim Documentation: https://docs.manim.community/en/stable/index.html
Superposition of Waves: https://github.com/RockSolid1106/Manim-Projects/tree/main/Superposition%20of%20waves
Pendulum: https://github.com/Uasoni/manim-pendulums
Transition in Hydrogen Wave Function: https://youtu.be/o3vFLcg10N0?si=OANL1alIhv6oftHf
Young's Double Slit Animation: https://youtu.be/sK1mU5omDoE?si=sh1TFsRpWRbociiW
Electron Double Slit Experiment: https://youtu.be/ZqS8Jjkk1HI?si=ex09lp7TjdSSaKFY
Ferminons Animation: https://github.com/quantum-visualizations/qmsolve
Gaussian Plane Plot: https://github.com/functionalvalue/manimexample.github.io/blob/main/07ParametricCurve04.py
Schrodinger Equation Simulation: https://youtu.be/v0UIGl4cTD0?si=GIZe-NsRuFdoQowe
Schrodinger Equation Solution for Hydrogen Atom: https://youtu.be/acN7E7AUHPk?si=9yXUrVH8WEfkonnp
Spot in Multiple Dimension, from Spiderman Across the Spiderverse: https://youtu.be/NO7_sRGSCng?si=4yQisVIpL-lzSUf_

"""

from manim import *
from manim.utils.rate_functions import ease_out_sine, ease_out_bounce
from math import sin, cos, sqrt, floor

import numpy as np
from numpy import square

# Intro Slide with names and such
class S1(Scene):
    def construct(self):
        intro = Text("FM132 - FOPW: Video Submission").scale(0.5).shift(UP*2).set_color(WHITE)
        self.play(FadeIn(intro))

        cop = Text("Copenhagen Interpretation &").scale(0.9).shift(UP*0.5).set_color_by_gradient(BLUE_B, LIGHT_BROWN)
        mws = Text("Many Worlds Theory").scale(0.9).next_to(cop, DOWN * 1.2).set_color_by_gradient(LIGHT_BROWN, BLUE_B)

        self.play(Write(VGroup(cop, mws)))

        names = VGroup(
            Text("Suraj Dayma - U20230102"),
            Text("Nikunj Agarwal - U20230068")
        ).arrange(DOWN).shift(DOWN*2).scale(0.3).set_color(BLUE_C)

        self.play(FadeIn(names))

        self.wait(1)

        self.play(FadeOut(VGroup(*[name for name in names], cop, mws, intro)))

        self.wait(1)

# Pendulum
class S2(Scene):
    def construct(self):
        time = ValueTracker(0)
        initThetaS = PI/4
        initTheta = Variable(initThetaS, r'\theta')
        l = 5
        g = 9.8
        w = np.sqrt(g/l)
        p = 2*PI/w

        originX = -2.25
        originY = 3
        startPoint = Dot([originX, originY, 0], radius=0.1, color=BLUE_B)
        originShft = originX* RIGHT + originY * UP

        theta = DecimalNumber().shift(10*UP)
        theta.add_updater(lambda m: m.set_value(initTheta.tracker.get_value()*np.sin(w*time.get_value())))
        self.add(theta)

        def getLine(x, y):
            line = Line(start=ORIGIN + originShft, end=x*RIGHT+y*UP+originShft, color=BLUE_B)
            global verticalLine
            verticalLine = DashedLine(start=line.get_start(), end=line.get_start() + l*DOWN, color=WHITE)

            return line

        def getEndBall(x, y):
            endBall = Dot(fill_color=BLUE_B, fill_opacity=1).move_to(x*RIGHT+y*UP+originShft).scale(l)
            return endBall

        line = always_redraw(lambda: getLine(l * np.sin(theta.get_value()), -l * np.cos(theta.get_value())))
        ball = always_redraw(lambda: getEndBall(l * np.sin(theta.get_value()), -l * np.cos(theta.get_value())))


        wDef = MathTex(r'w = \sqrt{\frac{g}{l}}').shift(UP*1, RIGHT*3)
        naturalFreqLabel = Text('"Natural Frequency"').scale(0.7).set_color_by_gradient(RED_B, YELLOW_B).align_to(wDef, LEFT).next_to(wDef, UP, buff=DEFAULT_MOBJECT_TO_MOBJECT_BUFFER+0.1).shift(1*RIGHT)
        wDef1 = MathTex(r'w = \sqrt{\frac{g}{l}}').shift(UP*1, RIGHT*3)
        thetaEquationMotion = MathTex(r'{{\theta}}(t)={{\theta}} \cos(wt)').set_color_by_tex("heta", BLUE).next_to(wDef, DOWN, buff=DEFAULT_MOBJECT_TO_MOBJECT_BUFFER+0.2).align_to(wDef, LEFT)
        thetaRange = MathTex(r'-3 < {{\theta}}< 3').set_color_by_tex("heta", BLUE).next_to(thetaEquationMotion, DOWN, buff=DEFAULT_MOBJECT_TO_MOBJECT_BUFFER+0.5).align_to(thetaEquationMotion, LEFT)


        pendulum = VGroup(startPoint, line, ball, verticalLine)
        text = VGroup(wDef1, naturalFreqLabel, thetaEquationMotion, wDef)


        self.play(FadeIn(pendulum))
        self.play(
            Write(text),
            time.animate.set_value(2*p), 
            rate_func=ease_out_sine, 
            run_time=5
            )

        self.play(VGroup(*pendulum, *text).animate.shift(DOWN*10))
        self.remove(*pendulum)

        self.wait()

# Quantum mechenics ka graph daaldo from top (S3) (if possible)
# No need to do this, video used.
class S3(Scene):
    def construct(self):
        pass

# Superposition of Waves
class S4(Scene):
    def construct(self):
        a = Axes(
            x_range=[-10, 10, PI/4],
            y_range=[-3, 3, PI/4], 
            y_length=3,
            axis_config={"include_tip":False, "include_ticks":False},
            color=WHITE
        ).move_to(UP*1.5)
        sin_g = a.plot(lambda x: np.sin(x), x_range=[-10*PI, -PI], color="#0E94ED") # Left wave
        cos_g = a.plot(lambda x: 0.5*np.sin(2*x), x_range=[PI, 10*PI], color="#ED2726") # Right wave
        
        
        
        tracker = ValueTracker(0)
        
        def get_sin_value(x): # Function to return the value of the left wave at an x coord. Returns 0 if input is outside the x_range of the wave.
            if x>=((-10*PI)+(PI*tracker.get_value())) and x<=(-PI+(PI*tracker.get_value())):
                return np.sin(x - PI*tracker.get_value())
            else:
                return 0
                
        def get_cos_value(x):
            if x>=(PI - (PI*tracker.get_value())) and x<=(10*PI - (PI*tracker.get_value())):
                return np.sin(2*x + 2*PI*tracker.get_value())/2
            else:
                return 0
    
        a2 = Axes(
            axis_config={"include_tip": False, "include_ticks":False},
            x_range=[0, 15],
            y_range=[-3, 3], 
            y_length=3,
        ).move_to(DOWN*1.5)
                
        tracer = Dot(a2.c2p(0, 0), color=GREEN)
        
        trace = TracedPath(tracer.get_center, stroke_color="#ED2726", stroke_width=5)

        self.play(Create(VGroup(a, a2, trace, tracer)))
        
        combination_g = a.plot(lambda x: get_sin_value(x) + get_cos_value(x), color="#0E94ED") # This is the visible graph in the upper graph. Its just the sum of values of the waves at all points.
        self.play(Create(combination_g))
        combination_g.add_updater(
            lambda x: x.become(a.plot(lambda x: get_sin_value(x) + get_cos_value(x), color="#0E94ED"))
        )
        
        tracer.add_updater(
            lambda x: x.move_to(a2.c2p(tracker.get_value(), get_sin_value(0) + get_cos_value(0)))
        )
        
        
        self.add(tracker)
        tracker.add_updater(lambda x, dt: x.increment_value(dt))
        
        self.wait(12)
        
        
        tracker.clear_updaters()

# Youngs Double Slit Animation - Didn't work out, just used the blender video
class S5(Scene):
    def construct(self):        

        black_rectangle = Rectangle(color=BLACK, fill_color=BLACK, fill_opacity=1).scale(10) #This Rectangle will fill into the entire screen

        a = ValueTracker(-2 * PI) #value of "a"
        axis = Axes()
        sin = axis.plot(lambda x : np.sin(x)).set_color(RED) #f(x) = sin(x)
        sin2 = always_redraw(lambda : axis.plot(lambda x : np.sin(x + a.get_value())).set_color(BLUE)) #f(x) = sin(x + a)
        sin3 = always_redraw(lambda : axis.plot(lambda x : np.sin(x + a.get_value()) + np.sin(x)).set_color(GREEN)) #f(x) = sin(x + a + sin(x))
        sin_label = Tex("$f(x) = sin(x)$").set_color(RED).to_corner(UL) #label for "f(x) = sin(x)"
        sin2_label = Tex("$f(x) = sin(x +  $").set_color(BLUE).to_edge(UP) #label for "f(x) = sin(x + a)"
        sin3_label = Tex("$f(x) = sin(x + sin(x) + $").set_color(GREEN).to_edge(DOWN).shift(LEFT) #label for f(x) = sin(x + a + sin(x))
        a_value = always_redraw(
            lambda : DecimalNumber(num_decimal_places=2).set_value(a.get_value()).set_color(BLUE).next_to(sin2_label, RIGHT)
        ) #Writes the value of "a" and is constantly updating
        a_copy = always_redraw(
            lambda : DecimalNumber(num_decimal_places=2).set_value(a.get_value()).set_color(GREEN).next_to(sin3_label, RIGHT)
        ) #Writes the value of "a" and is constantly updating
        sin2_brac = Tex("$)$").next_to(a_value, RIGHT).set_color(BLUE)
        sin3_brac = Tex("$)$").next_to(a_copy, RIGHT).set_color(GREEN)

        self.play(Create(axis))
        self.play(Write(VGroup(sin_label, sin2_label, sin3_label)))
        self.play(Write(VGroup(a_value, a_copy, sin2_brac, sin3_brac)))
        self.play(Create(VGroup(sin, sin2, sin3)))
        self.wait()
        self.play(a.animate.set_value(5 * PI), run_time=10) #Value of "a" becomes (5 * PI)
        self.wait()
        self.play(a.animate.set_value(-5 * PI), run_time=20) #Value of "a" becomes (-5 * PI)
        self.wait()
        self.play(FadeIn(black_rectangle)) #This animation will have an effect of every mobject on screen "Fading Out"
        self.wait()

# Light is a wave and particle both
class S6(Scene):
    def construct(self):
        light = Text("Light").scale(0.8).set_color(YELLOW)

        self.play(FadeIn(light))

        self.play(light.animate.shift(UP*2))

        wave_arrow = Arrow(
            light.get_bottom(), 
            light.get_bottom() + [-2, -2, 0], 
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        particle_arrow = Arrow(
            light.get_bottom(),
            light.get_bottom() + [2, -2, 0],
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        arrows = VGroup(wave_arrow, particle_arrow)

        self.play(DrawBorderThenFill(arrows))

        wave_text = Text("Wave").scale(0.7).next_to(wave_arrow.get_end(), DOWN).set_color(LIGHT_BROWN)
        particle_text = Text("Particle").scale(0.7).next_to(particle_arrow.get_end(), DOWN).set_color(LIGHT_BROWN)

        texts = VGroup(wave_text, particle_text)

        axes = Axes(
            x_range=[-10, 10, 1],
            y_range=[-10, 10, 10],
            x_length=2,
        )

        sin_graph = axes.plot(lambda x: np.sin(x), color=BLUE).rotate(-PI/2).next_to(wave_text, DOWN)

        cols = [BLUE_A, BLUE_B, BLUE_C, BLUE_D, BLUE_E, PURPLE_A, PURPLE_B]

        circle = Circle(radius=0.5, color=GREEN)
        num_points = 6
        angles = [n * (360 / num_points) for n in range(num_points)]
        points = [circle.point_at_angle(n*DEGREES) for n in angles]
        circles = VGroup(*[Dot(fill_opacity=0.6, color=cols[p]).move_to(points[p]) for p in range(0, len(points))]).next_to(particle_text, DOWN)

        self.play(FadeIn(texts))
        self.play(Create(sin_graph), Create(circles))

        self.wait()

        all_light = VGroup(*texts, sin_graph, *circles, *arrows, light)
        self.play(all_light.animate.shift(LEFT*3.5))

        divider_line = Line(UP*10, DOWN*10).set_color(GRAY)
        self.play(Create(divider_line))

        GPThompsonPic = ImageMobject("George_Paget_Thomson.jpg").scale(1).next_to(divider_line, RIGHT, buff=.7)
        DavGermPic = ImageMobject("Davisson_and_Germer.jpg").scale(1.15).next_to(GPThompsonPic, RIGHT, buff=.8)
        self.play(FadeIn(GPThompsonPic), FadeIn(DavGermPic))

        GPText = Text("GP Thompson").scale(0.3).next_to(GPThompsonPic, DOWN)
        DGText = Text("Davisson & Germe").scale(0.3).next_to(DavGermPic, DOWN)

        self.play(FadeIn(VGroup(GPText, DGText)))

        self.wait()

        self.play(
            Group(*all_light, divider_line).animate.shift(LEFT * 10), 
            GPThompsonPic.animate.shift(LEFT * 7), 
            DavGermPic.animate.shift(LEFT * 7),
            VGroup(GPText, DGText).animate.shift(LEFT * 7)
            )
        
        electron = LabeledDot(MathTex(r"e^-"), fill_opacity = 0.3, stroke_width=1, color=RED).shift(UP * 2).shift(RIGHT * 3)

        self.play(FadeIn(electron))

        wave_arrow = Arrow(
            electron.get_bottom(), 
            electron.get_bottom() + [-2, -2, 0], 
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        particle_arrow = Arrow(
            electron.get_bottom(),
            electron.get_bottom() + [2, -2, 0],
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        arrows = VGroup(wave_arrow, particle_arrow)

        self.play(DrawBorderThenFill(arrows))

        wave_text = Text("Wave").scale(0.7).next_to(wave_arrow.get_end(), DOWN).set_color(LIGHT_BROWN)
        particle_text = Text("Particle").scale(0.7).next_to(particle_arrow.get_end(), DOWN).set_color(LIGHT_BROWN)

        texts = VGroup(wave_text, particle_text)

        axes = Axes(
            x_range=[-10, 10, 1],
            y_range=[-10, 10, 10],
            x_length=2,
        )

        sin_graph = axes.plot(lambda x: np.sin(x), color=BLUE).rotate(-PI/2).next_to(wave_text, DOWN)

        cols = [BLUE_A, BLUE_B, BLUE_C, BLUE_D, BLUE_E, PURPLE_A, PURPLE_B]

        circle = Circle(radius=0.5, color=GREEN)
        num_points = 6
        angles = [n * (360 / num_points) for n in range(num_points)]
        points = [circle.point_at_angle(n*DEGREES) for n in angles]
        circles = VGroup(*[Dot(fill_opacity=0.6, color=cols[p]).move_to(points[p]) for p in range(0, len(points))]).next_to(particle_text, DOWN)

        self.play(FadeIn(texts))
        self.play(Create(sin_graph), Create(circles))

        self.wait(2)

        matter = Text("Matter").scale(0.8).set_color(PINK).move_to(electron.get_center())
        self.play(ReplacementTransform(electron, matter))

        self.wait(2)

# Citation of Electron Double Slit Experiment
class S_Citation(Scene):
    def construct(self):
        bg = Rectangle(color=GREEN, fill_opacity = 1, height=20, width=20)
        self.add(bg)

        citation = Text("Source: (Bach et al., Controlled double-slit electron diffraction 2013)", slant=ITALIC).scale(0.4).to_corner(DL)
        citation.set_color(YELLOW_B)
        self.add(citation)

# Schrodinger Cat Experiment Text
class S7_1(Scene):
    def construct(self):
        text = Text("Schrödinger's Cat", slant=ITALIC).set_color_by_gradient(RED_D, LIGHT_BROWN).scale(1.3)
        self.play(Write(text), run_time=1.5)
        self.wait()
        self.play(FadeOut(text))
        self.wait()

# Schrodinger Cat Setup
class S7(Scene):
    def construct(self):
        catAlive = SVGMobject("catAlive.svg")
        catDead = SVGMobject("catDead.svg")

        one = SVGMobject("1.svg")
        three = SVGMobject("3.svg")
        two = SVGMobject("2.svg")
        four = SVGMobject("4.svg")
        five = SVGMobject("5.svg")
        six = SVGMobject("6.svg")

        self.play(DrawBorderThenFill(one))
        self.play(FadeIn(two))
        self.play(FadeIn(three))
        self.play(FadeIn(four))
        self.play(FadeIn(five))
        self.play(FadeIn(six))

        self.wait(2)

        self.remove(one, two, three, four, five, six)
        self.add(catAlive)

        self.play(
            catAlive.animate.shift(LEFT*3)
        )

        arrow_dead = Arrow(
            catAlive.get_right(),
            catAlive.get_right() + [2, 2, 0],
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        arrow_alive = Arrow(
            catAlive.get_right(),
            catAlive.get_right() + [2, -2, 0],
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        catDead.next_to(arrow_dead.get_end(), RIGHT)
        catAlive2 = catAlive.copy()
        catAlive2.next_to(arrow_alive.get_end(), RIGHT)

        self.play(
            Create(arrow_dead)
        )

        self.play(
            FadeIn(catDead)
        )

        self.play(
            Create(arrow_alive)
        )
        
        self.play(
            FadeIn(catAlive2)
        )

        self.wait(2)
        
        self.play(
            FadeOut(*self.mobjects)
        )

        self.wait()

# Atom is constantly in superposition
class S8(Scene):
    def construct(self):
        notDecayed = Circle(0.2, color=BLUE, stroke_color=GRAY, stroke_width=0.5, fill_opacity=0.3).scale(4).shift(LEFT*2)
        notDecayedLabel = Text("Not Decayed", color=BLUE).scale(0.4).next_to(notDecayed, DOWN)
        safe = VGroup(notDecayed, notDecayedLabel)
        notDecayedBox = SurroundingRectangle(notDecayed, stroke_width=0.8)
        
        decayed = Circle(0.2, color=RED, stroke_color=GRAY, stroke_width=0.5, fill_opacity=0.3).scale(4).shift(RIGHT*2)
        decayedLabel = Text("Decayed", color=RED).scale(0.4).next_to(decayed, DOWN)
        risky = VGroup(decayed, decayedLabel)
        decayedBox = SurroundingRectangle(decayed, stroke_width=0.8)

        self.play(
            Create(notDecayedBox),     
            Create(decayedBox),
        )

        self.play(
            DrawBorderThenFill(notDecayed),
            Create(notDecayedLabel),
            DrawBorderThenFill(decayed),
            Create(decayedLabel),
        )


        text2 = Text("State?", color=BLUE).scale(0.4).shift(DOWN * 1.15)
        text3 = Text("", color=BLUE).scale(0.4).shift(DOWN * 1.20)

        self.play(
            decayed.animate.move_to(ORIGIN),
            notDecayed.animate.move_to(ORIGIN),
            decayedBox.animate.move_to(ORIGIN),
            notDecayedBox.animate.move_to(ORIGIN),
            Transform(notDecayedLabel, text2),
            Transform(decayedLabel, text3),
        )

        self.wait(2)

        self.play(VGroup(safe, risky, notDecayedBox, decayedBox).animate.shift(LEFT * 4))

        heading = Text("Schrödinger Equation", color=LIGHT_BROWN).scale(0.7).shift((RIGHT + UP) * 2.5)
        schrodinger_eq = MathTex(r"i \hbar \frac{\partial}{\partial t}\Psi(\mathbf{r},t) = \hat H \Psi(\mathbf{r},t)", color=YELLOW).scale(0.7).next_to(heading, DOWN)
        self.play(FadeIn(heading))
        self.play(Write(schrodinger_eq))

        self.wait(2)

        self.play(
            VGroup(safe, risky, notDecayedBox, decayedBox).animate.shift(DOWN * 1)
        )

        arrow1 = Arrow(decayedBox.get_right(), decayedBox.get_right() + ((RIGHT * 2) + (UP * 1)), 
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,)
        arrow2 = Arrow(decayedBox.get_right(), decayedBox.get_right() + ((RIGHT*2 + DOWN*1)), 
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,)
        
        blueDot = Dot(color=BLUE).scale(3).next_to(arrow1.get_end(), RIGHT)
        blueText = Text("Not Decayed.").scale(0.4).next_to(blueDot, RIGHT)

        redDot = Dot(color=RED).scale(3).next_to(arrow2.get_end(), RIGHT)
        redText = Text("Decayed.").scale(0.4).next_to(redDot, RIGHT)

        s1 = notDecayedLabel.copy()
        s2 = s1.copy()

        self.play(
            Create(arrow1),
            Create(arrow2)
        )

        self.add(s1,s2)
        self.play(
            DrawBorderThenFill(blueDot),
            DrawBorderThenFill(redDot),
            ReplacementTransform(s1, blueText),
            ReplacementTransform(s2, redText)
        )

        self.wait(2)

# While the box is closed, cat is both dead and alive
class S9(Scene):
    def construct(self):
        catAlive = SVGMobject("catAlive.svg", fill_color=BLUE, fill_opacity=0.1).shift(LEFT * 2)
        catDead = SVGMobject("catDead.svg", fill_color=RED, fill_opacity=0.1).shift(RIGHT * 2)

        alive = Text("Alive", color=BLUE).scale(0.4).next_to(catAlive, DOWN)
        dead = Text("Dead", color=BLUE).scale(0.4).next_to(catDead, DOWN)

        s1 = Text("Alive AND Dead", color=RED).scale(0.4).move_to([0, catAlive.get_bottom()[1] - 0.3, 0])
        s2 = Text("")

        self.play(
            FadeIn(catAlive),
            FadeIn(catDead),
            FadeIn(alive),
            FadeIn(dead)
            )
        
        self.wait(2)

        self.play(
            catAlive.animate.move_to(ORIGIN),
            catDead.animate.move_to(ORIGIN),
            Transform(alive, s1),
            Transform(dead, s2)
        )

        self.wait(2)

        self.play(FadeOut(*self.mobjects))

        self.wait(1)

# Wave function 
class S10(Scene):
    def construct(self):

        heading = Text("Wave Function").to_edge(UP, buff=1.3).scale(0.8).set_color(YELLOW_A)
        self.play(FadeIn(heading))

        wave_func = MathTex(r"\Psi(x, t)", color=WHITE)
        wave_func[0][0:1].set_color(YELLOW)
        self.play(Write(wave_func))

        self.wait()

        self.play(wave_func.animate.shift(UP * 0.6))

        wave_sq = MathTex(r"|{\Psi(x, t)|^2").next_to(wave_func, DOWN)
        wave_sq[0][1:2].set_color(YELLOW)
        self.play(Write(wave_sq))

        self.wait()

        self.play(wave_sq.animate.shift(LEFT*0.8))
        prob = MathTex(r"= \rho(x)").next_to(wave_sq, RIGHT,)
        prob[0][1:2].set_color(GREEN_C)
        self.play(Write(prob))

        blueDot = Circle(0.6, color=BLUE, stroke_color=GRAY, stroke_width=0.5, fill_opacity=0.3).next_to(prob, DOWN*3 + LEFT*7)
        redDot = Circle(0.6, color=RED, stroke_color=GRAY, stroke_width=0.5, fill_opacity=0.3).next_to(blueDot, RIGHT*2)


        self.play(
            DrawBorderThenFill(blueDot),
            DrawBorderThenFill(redDot),
        )

        self.play(
            redDot.animate.move_to(blueDot.get_center())
        )

        self.wait()

        catAlive = SVGMobject("catAlive.svg").next_to(redDot, RIGHT, buff = 1)
        catDead = SVGMobject("catDead.svg").next_to(redDot, RIGHT, buff = 1)

        alive = Text("Alive", color=BLUE).scale(0.4).next_to(catAlive, DOWN)
        dead = Text("Dead", color=BLUE).scale(0.4).next_to(catDead, DOWN)

        self.play(
            FadeIn(catAlive),
            FadeIn(alive)
        )

        self.wait()

        self.play(
            ReplacementTransform(catAlive, catDead),
            ReplacementTransform(alive, dead)
        )

        self.wait(2)

        self.play(FadeOut(*self.mobjects))

        self.wait(2)

# Cop. Inter. Text Green Screen
class S12(Scene):
    def construct(self):
        bg = Rectangle(color=GREEN, fill_opacity = 1, height=20, width=20)
        self.add(bg)

        text = Text("Copenhagen Interpretation", slant=ITALIC)
        text.set_color_by_gradient(BLUE_B, LIGHT_BROWN)
        self.play(Write(text))
        self.wait(8)


# Trying to Generate 3D Gaussian Surface which would transform into 
# a single point upon "measurement" -> didn't work
class Gaussian(ThreeDScene):
    def construct(self):
        axis_config = {
            "x_range": [-5, 5],
            "y_range": [-5, 5],
            "z_range": [-5, 5],
        }

        axis = ThreeDAxes()
        self.set_camera_orientation(phi=50 * DEGREES, theta=60 * DEGREES, focal_distance=20)
        self.add(axis)
        self.begin_ambient_camera_rotation(rate=1)

        u0 = 0
        v0 = 0
        sig = 0.8
        res = (30, 30)

        curve = Surface(
            lambda u, v:
            np.array([
                u, v,
                2.5 * np.exp(
                    -0.1 * (( (u - u0)/sig )**2 + ( (v - v0) / sig )**2)
                )
            ]),
            u_range=[-TAU + 3, TAU - 3],
            v_range=[-TAU + 3, TAU - 3],
            resolution=res,
            checkerboard_colors=[BLUE, PINK]
        )

        sig2 = 0.01
        curve2 = Surface(
            lambda u, v:
            np.array([
                u, v,
                2.5 * np.exp(
                    -0.1 * (( (u - u0)/sig2 )**2 + ( (v - v0) / sig2 )**2)
                )
            ]),
            u_range=[-TAU + 3, TAU - 3],
            v_range=[-TAU + 3, TAU - 3],
            resolution=res,
            checkerboard_colors=[PINK, PURPLE]
        )

        self.play(Write(curve), run_time=4)
        self.wait(7)
        self.play(ReplacementTransform(curve, curve2))
        self.wait(6)


# Start of Many Worlds and Entanglement
class M1(Scene):
    def construct(self):
        catAlive = SVGMobject("catAlive.svg", fill_color=BLUE, fill_opacity=0.1).shift(LEFT * 2)
        catDead = SVGMobject("catDead.svg", fill_color=RED, fill_opacity=0.1).shift(RIGHT * 2)

        alive = Text("Alive", color=BLUE).scale(0.4).next_to(catAlive, DOWN)
        dead = Text("Dead", color=BLUE).scale(0.4).next_to(catDead, DOWN)

        s1 = Text("Alive AND Dead", color=RED).scale(0.4).move_to([0, catAlive.get_bottom()[1] - 0.3, 0])
        s2 = Text("")

        self.play(
            FadeIn(catAlive),
            FadeIn(catDead),
            FadeIn(alive),
            FadeIn(dead)
            )
        

        self.wait(2)

        aliveGroup = VGroup(catAlive, alive)
        deadGroup = VGroup(catDead, dead)

        self.play(
            aliveGroup.animate.shift(UP * 2),
            deadGroup.animate.shift(UP * 2)
        )

        notDecayed = Circle(0.2, color=BLUE, stroke_color=GRAY, stroke_width=0.5, fill_opacity=0.3).scale(4).shift(LEFT*2)
        notDecayedLabel = Text("Not Decayed", color=BLUE).scale(0.4).next_to(notDecayed, DOWN)
        safe = VGroup(notDecayed, notDecayedLabel).scale(0.7).shift(DOWN*2)
        
        decayed = Circle(0.2, color=RED, stroke_color=GRAY, stroke_width=0.5, fill_opacity=0.3).scale(4).shift(RIGHT*2)
        decayedLabel = Text("Decayed", color=RED).scale(0.4).next_to(decayed, DOWN)
        risky = VGroup(decayed, decayedLabel).scale(0.7).shift(DOWN*2)

        a1 = Arrow(
            aliveGroup.get_bottom(), 
            safe.get_top(),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )
        a2 = Arrow(
            deadGroup.get_bottom(), 
            risky.get_top(),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        self.play(
            DrawBorderThenFill(a1),
            DrawBorderThenFill(a2),
            DrawBorderThenFill(safe),
            DrawBorderThenFill(risky)
        )

        self.wait(3)

        text2 = Text("State?", color=BLUE).scale(0.4).shift(DOWN * 1.15)
        text3 = Text("", color=BLUE).scale(0.4).shift(DOWN * 1.20)

        self.play(
            FadeOut(aliveGroup, deadGroup, a1, a2),
        )

        self.play(
            decayed.animate.move_to(ORIGIN),
            notDecayed.animate.move_to(ORIGIN),
            Transform(notDecayedLabel, text2),
            Transform(decayedLabel, text3),
        )

        atom = VGroup(decayed, notDecayed, notDecayedLabel)
        self.play(
            atom.animate.shift(LEFT * 5)
        )

        detector = SVGMobject("detector.svg", fill_color=YELLOW, fill_opacity=0.6).scale(0.6).shift(LEFT * 2)

        a1 = Arrow(
            notDecayed.get_right(), 
            detector.get_left(),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        self.play(
            DrawBorderThenFill(detector),
            DrawBorderThenFill(a1)
        )
        
        hammer = SVGMobject("4.svg", fill_color=ORANGE, fill_opacity=0.6).scale(0.6).shift(RIGHT * 1)

        a2 = Arrow(
            detector.get_right(), 
            hammer.get_left(),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        self.play(
            DrawBorderThenFill(hammer),
            DrawBorderThenFill(a2)
        )

        cat = SVGMobject("2.svg", fill_color=GRAY, fill_opacity=0.6).scale(0.6).shift(RIGHT * 4.5)

        a3 = Arrow(
            hammer.get_right(), 
            cat.get_left(),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        self.play(
            DrawBorderThenFill(cat),
            DrawBorderThenFill(a3)
        )

        self.wait(2)

        text = Text("Quantum Entanglement", slant=ITALIC).shift(UP*2)
        text.set_color_by_gradient(BLUE_B, LIGHT_BROWN)
        self.play(Write(text))
        self.wait(5)
        
        self.play(FadeOut(atom, a1, detector, a2, hammer, a3, cat))

        # Quantum Entanglement 

        point1 = Text("Occurs at Atomic Level")
        point2 = Text("Singular Wave Function for the World")
        
        points = VGroup(point1, point2)


        # points aesthetics
        points.scale(0.5).set_color(YELLOW, YELLOW_D)

        # points position 
        points.arrange(DOWN, buff=0.3) #, center=False, aligned_edge=LEFT)

        for point in points:
            self.play(FadeIn(point))
            self.wait()
        
        self.wait()
        self.play(FadeOut(*self.mobjects))

        self.wait(1)


        detector = SVGMobject("detector.svg", fill_color=YELLOW, fill_opacity=0.6).scale(0.6).shift(LEFT * 2)


        self.play(
            DrawBorderThenFill(atom),
            DrawBorderThenFill(detector),
            DrawBorderThenFill(a1),
            DrawBorderThenFill(hammer),
            DrawBorderThenFill(a2),
            DrawBorderThenFill(cat),
            DrawBorderThenFill(a3),
        )

        self.wait(2)

        self.play(
            atom.animate.shift(LEFT * 15),
            detector.animate.shift(LEFT * 15),
            a1.animate.shift(LEFT * 15),
            hammer.animate.shift(LEFT * 15),
            a2.animate.shift(LEFT * 15),
            cat.animate.shift(LEFT * 9.5),
            a3.animate.shift(LEFT * 15),
        )

        self.wait()

        air = SVGMobject("air.svg", fill_color=BLUE_A, fill_opacity=0.6).scale(0.5).next_to(cat, buff = RIGHT * 2)
        a1 = Arrow(
            cat.get_right(), 
            air.get_left(),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        self.play(
            DrawBorderThenFill(air),
            DrawBorderThenFill(a1)
        )

        box = SVGMobject("box.svg", fill_color=GRAY, fill_opacity=0.3).scale(0.6).next_to(air, buff = RIGHT * 2)
        a2= Arrow(
            air.get_right(), 
            box.get_left(),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        self.play(
            DrawBorderThenFill(box),
            DrawBorderThenFill(a2)
        )

        self.wait(2)

        you = Text("You", slant=ITALIC, color=WHITE).scale(0.7).next_to(box, buff = RIGHT * 2)
        a3= Arrow(
            box.get_right(), 
            you.get_left(),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        self.play(
            DrawBorderThenFill(you),
            DrawBorderThenFill(a3)
        )

        self.wait(2)

        self.play(
            FadeOut(cat, air, box, a1, a2, a3),
            you.animate.move_to(ORIGIN + UP),
        )

        decayed.move_to(ORIGIN + DOWN)
        notDecayed.move_to(ORIGIN + DOWN)
        self.play(DrawBorderThenFill(decayed), DrawBorderThenFill(notDecayed))

        a1= Arrow(
            you.get_bottom(), 
            notDecayed.get_top() + np.array([-1, 0, 0]),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )
        a2= Arrow(
            you.get_bottom(), 
            decayed.get_top() + np.array([1, 0, 0]),
            max_tip_length_to_length_ratio=0.1,
            stroke_width=3,
            color=YELLOW_B,
        )

        self.play(
            notDecayed.animate.shift(LEFT),
            decayed.animate.shift(RIGHT),
            DrawBorderThenFill(a1),
            DrawBorderThenFill(a2),
        )

        self.wait(2)

        self.play(
            FadeOut(*self.mobjects)
        )


# Quantum Entanglement Text
class M1_Text(Scene):
    def construct(self):
        bg = Rectangle(color=GREEN, fill_opacity = 1, height=20, width=20)
        self.add(bg)

        text = Text("Quantum Entanglement", slant=ITALIC)
        text.set_color_by_gradient(BLUE_B, LIGHT_BROWN)
        self.play(Write(text))
        self.wait(5)

# There is another you, many worlds animation
# Scene starts with a person and then zooms out 
# to reveal multiple verisons
class M2(ZoomedScene):
    def construct(self):
        people = VGroup()
        for i in range(200):
            person = SVGMobject("person.svg", fill_color=random_bright_color(), fill_opacity=0.6).scale(0.4)
            people.add(person)

        people.arrange_in_grid(10, 20, buff=0.7)

        self.camera.frame.save_state()

        start = 109
        end = start + 2

        initial = VGroup(people[start:end])

        self.play(self.camera.frame.animate.set(
            width=initial.width + 2,
            height=initial.height + 2,
            ).move_to(initial))
        self.play(
            DrawBorderThenFill(initial),
        )
        self.wait(2)
        
        people = VGroup(*[people[i] for i in range(len(people)) if i not in [start, start + 1]])

        self.play(
            self.camera.frame.animate.set(width=people.width + 0.3, height=people.height + 0.3).move_to(people.get_center()),
            FadeIn(people, run_time=5),
        )

        self.wait(2)

# Many Worlds Interpretation Text
class M2_Text(Scene):
    def construct(self):
        bg = Rectangle(color=GREEN, fill_opacity = 1, height=20, width=20)
        self.add(bg)

        text = Text("Many Worlds Interpretation", slant=ITALIC)
        text.set_color_by_gradient(BLUE_B, LIGHT_BROWN)
        self.play(Write(text))
        self.wait(5)
        
class Outro(Scene):
    def construct(self):
        pic = ImageMobject("images/Erwin_Schrödinger_(1933).jpg").scale(1.4).shift(LEFT * 4 + UP * 0.5)
        
        quote = Paragraph("\"I don't like it, and I'm sorry \nI ever had anything to do with it.\"", alignment='center')
        quote.scale(0.7).next_to(pic, RIGHT * 3).set_color_by_gradient(BLUE_B, LIGHT_BROWN)

        name = Text(" - Erwin Schrödinger", color=YELLOW).scale(0.4).next_to(quote, DOWN)         
        on = Paragraph("on probability interpretation\n of quantum mechanics",  alignment='center', color=YELLOW_D)
        on.scale(0.3).next_to(name, DOWN * 0.6)

        self.play(
            FadeIn(pic),
            Write(quote)
            )
        
        self.play(
            FadeIn(name),
            FadeIn(on)
        )
        

        self.wait()

