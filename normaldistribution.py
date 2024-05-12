# Reference:
# https://github.com/thomasnield/3mds/blob/main/20230822_manim_normal_distribution_walkthrough.py#L30

from manim import *
from scipy.stats import norm
import os

data = np.array([65.27153711, 61.69996242, 60.98565375, 65.30031155, 63.51806848, 68.19351011
                    , 66.95478689, 64.55759847, 63.39196506, 67.54289154, 63.19717054, 67.49928145
                    , 63.19766386, 72.39460819, 65.06618895, 61.47292356, 63.65363793, 67.40224834
                    , 69.5453564, 60.68221867, 61.8700392, 61.83843595, 68.61144043, 65.56311344
                    , 64.75033447, 58.11010178, 58.05563835, 64.47898656, 61.58413871, 62.73034603
                    , 71.43664684, 63.16241713, 64.67237221, 57.72834468, 69.02743165, 61.51629892
                    , 63.11713525, 66.89447031, 64.641517, 64.27097855, 68.33858714, 61.08831661
                    , 70.21508249, 61.28410416, 66.28841856, 66.36769216, 66.85034888, 60.67420648
                    , 60.26603058, 65.16350301, 65.50142234, 63.43139858, 65.77102549, 62.59483884
                    , 64.99651396, 63.70706524, 63.53666649, 63.56936279, 67.49544015, 65.61582933
                    , 70.70856991, 64.50069849, 58.86870671, 66.22804048, 65.12576775, 58.25042313
                    , 60.85683308, 63.962537, 63.76654245, 62.21488775, 66.89233257, 66.27436643
                    , 66.06183922, 62.58999392, 62.01327474, 69.776555, 65.86363553, 66.37112032
                    , 65.0016078, 68.24377827, 60.02304379, 64.70747144, 62.58229384, 69.80994781
                    , 67.92226975, 62.4487201, 63.18153599, 65.53453952, 65.39880782, 59.37181606
                    , 67.58819312, 65.053418, 62.32011733, 65.51449093, 61.70972692, 66.08806211
                    , 63.49776397, 68.8884009, 63.55453324, 66.02846214])


class PDFtoCDFtoPPFScene(MovingCameraScene):

    def construct(self):
        skip_animations = False

        # Declare PDF model
        class PDFPlot(VGroup):

            def __init__(self, mean, std):
                super().__init__()

                f = lambda x: norm.pdf(x, mean, std)

                self.mean = mean
                self.std = std
                self.lower_x = mean - std * 3
                self.upper_x = mean + std * 3

                axes = Axes(x_range=[self.lower_x, self.upper_x, std],
                            y_range=[0, f(mean) + .1, (f(mean) + .1) / 4],
                            x_axis_config={"include_numbers": False,
                                           "numbers_to_exclude": [mean - 4 * std]
                                           },
                            y_axis_config={"include_numbers": True,
                                           "decimal_number_config": {
                                               "num_decimal_places": 2
                                           }
                                           }
                            )

                plot = axes.plot(f, color=BLUE)
                self.add(axes, plot)

                self.f = f
                self.axes = axes
                self.plot = plot

            def x2p(self, x):
                return self.axes.c2p(x, self.f(x))

            def area_range(self, x_start, x_end, color=BLUE):
                return self.axes.get_area(self.plot, color=color, x_range=(x_start, x_end))

        # Declare CDF model
        class CDFPlot(VGroup):
            def __init__(self, mean, std):
                super().__init__()
                f = lambda x: norm.cdf(x, mean, std)

                axes = Axes(x_range=[mean - 3 * std, mean + 3 * std, std],
                            y_range=[0, 1.1, .25],
                            y_axis_config={"include_numbers": True,
                                           "decimal_number_config": {
                                               "num_decimal_places": 2
                                           }
                                           }
                            )

                plot = axes.plot(f, color=RED)

                self.add(axes, plot)
                self.f = f
                self.axes = axes
                self.plot = plot

            def x2p(self, x):
                return self.axes.c2p(x, self.f(x))

            def plot_to_x(self, x):
                return self.axes.plot(self.f, color=RED, x_range=[mean - 3 * std, x])

            def vertical_line(self, x):
                return DashedLine(
                    start=self.axes.c2p(x, self.f(x)),
                    end=self.axes.c2p(x, 0),
                    color=RED
                )

            def horizontal_line(self, x):
                return DashedLine(
                    start=self.axes.c2p(-3 * std + mean, self.f(x)),
                    end=self.axes.c2p(x, self.f(x)),
                    color=RED
                )

        # Declare PPF model
        class PPFPlot(VGroup):
            def __init__(self, mean, std):
                super().__init__()
                f = lambda x: norm.ppf(x, mean, std)

                axes = Axes(x_range=[.001, .999, .05],
                            y_range=[mean - 3 * std, mean + 3 * std, std],
                            x_length=4,
                            x_axis_config={
                                "numbers_to_include": [0, .25, .5, .75, 1]
                            }
                            )

                plot = axes.plot(f, color=ORANGE, use_smoothing=True)

                self.add(axes, plot)
                self.f = f
                self.axes = axes
                self.plot = plot

            def p2p(self, p):
                return self.axes.c2p(p, self.f(p))

            def vertical_line(self, p):
                return DashedLine(
                    start=self.axes.c2p(p, -3 * std + mean),
                    end=self.axes.c2p(p, self.f(p)),
                    color=RED
                )

            def horizontal_line(self, p):
                return DashedLine(
                    start=self.axes.c2p(p, self.f(p)),
                    end=self.axes.c2p(0, self.f(p)),
                    color=RED
                )

        mean, std = data.mean(), data.std()

        # create PDF and CDF
        self.next_section("Create PDF and CDF, and zoom out", skip_animations=skip_animations)
        pdf_model = PDFPlot(mean, std)
        cdf_model = CDFPlot(mean, std)

        # stack the PPF to the right for later
        ppf_model = PPFPlot(mean, std).to_edge(RIGHT)

        # stack the PDF and CDF to the left
        left_panel = VGroup(cdf_model, pdf_model) \
            .arrange(DOWN) \
            .scale_to_fit_height(7)

        # add PDF to the scene
        # self.add(pdf_model)

        # start camera zoomed on PDF
        self.camera.frame.save_state()
        self.camera.frame.scale(0.6).move_to(pdf_model)
        self.play(
            FadeIn(pdf_model, run_time=3)
        )
        self.wait()

        # zoom out the camera to reveal the CDF model axes
        self.play(
            Restore(self.camera.frame),
            FadeIn(cdf_model.axes)
        )
        self.wait()

        self.next_section("Draw area between PDF and CDF, label both", skip_animations=skip_animations)

        # Declare the range for x values on PDF
        x_upper_tracker = ValueTracker(mean - std * 3)

        # Declare the area for the PDF which will update based on the trackers above
        area_color = BLUE
        area: Mobject = always_redraw(
            lambda: pdf_model.area_range(-3 * std + mean, x_upper_tracker.get_value(), color=area_color)
        )

        # Draw the connecting dashed line between the PDF and CDF projecting the area
        connecting_line: DashedLine = always_redraw(lambda: DashedLine(
            start=pdf_model.x2p(x_upper_tracker.get_value()),
            end=cdf_model.x2p(x_upper_tracker.get_value()),
            color=RED
        ))

        # Project the area to the CDF as x_upper increases, also show the area as a decimal
        cdf_partial_plot = always_redraw(lambda: cdf_model.plot_to_x(x_upper_tracker.get_value()))
        area_label = always_redraw(lambda: DecimalNumber(cdf_model.f(x_upper_tracker.get_value()), num_decimal_places=2) \
                                   .scale(.8) \
                                   .next_to(cdf_model.x2p(x_upper_tracker.get_value()), RIGHT)
                                   )

        # Populate the area plot, connecting line, partial CDF plot, area label
        self.play(*[FadeIn(mobj) for mobj in (area, connecting_line, cdf_partial_plot, area_label)])

        # Run the animation for 5 seconds by filling the whole PDF/CDF
        self.play(
            x_upper_tracker.animate.set_value(mean + 3 * std),
            run_time=5
        )
        self.wait()
        self.wait()



# Execute rendering
if __name__ == "__main__":
    os.system(r"manim -qk -v WARNING -p --disable_caching -o PDFtoCDFtoPPFScene.mp4 082220223_manim_normal_distribution_walkthrough.py.py PDFtoCDFtoPPFScene")