module Utils

using FFTW
using Plots
using Images
using TestImages
using Interact

export fftplot, fftimg, fftimggui, seriesplot

"""
    seriesplot(x, s; kwargs...)

Plot a discrete signal `s` on coordinates `x`.
"""
function seriesplot(x, s; kwargs...)
    xwidth = (maximum(x) - minimum(x)) / length(x)
    p = bar(x, s; bar_width=0.1*xwidth, fillcolor=:black, legend=false, kwargs...)
    plot!(x, s; markershape=:circle, line=false, markerstrokewidth=0, legend=false)
    return(p)
end

"""
    fftplot(x, s)

Plot the signal `s` over `x`, as well as its amplitude and frequency spectrum.
"""
function fftplot(x, s)
    # plot the original signal
    p1 = seriesplot(x, s, title="signal")
    # calculate Fourier transform
    S = fftshift(fft(s))
    freqs = LinRange(-length(x)//2, length(x)//2, length(x))
    # plot amplitude and phase spectra
    p2 = seriesplot(freqs, abs.(S), title="amplitude")
    p3 = seriesplot(freqs, angle.(S), title="phase")
    # arrange plot layout
    l = @layout [a{0.4w} [b; c]]
    return(plot(p1, p2, p3, layout=l))
end

"""
    fftimg(img)
    
Display image `img` as well as its Fourier transform.
"""
function fftimg(img)
    p1 = heatmap(img, color=:grays, xticks=false, yticks=false, title="image")
    spectrum = fftshift(fft(channelview(img)))
    p2 = heatmap(log.(1 .+ abs.(spectrum)), color=:grays, xticks=false, yticks=false, title="amplitude")
    p3 = heatmap(angle.(spectrum), color=:grays, xticks=false, yticks=false, title="phase")
    l = @layout [a{0.4w} b{0.4h} c]
    return(plot(p1, p2, p3, layout=l))
end

function fftimshow(img; title="image")
    heatmap(img, color=:grays, xticks=false, yticks=false, title=title)
end
function fftpower(img; title="amplitude", logscale=true)
    if logscale
        spectrum = log.(1 .+ abs.(fftshift(fft(img))))
    else
        spectrum = abs.(fftshift(fft(img)))
    end
    heatmap(spectrum, color=:grays, xticks=false, yticks=false, title=title)
end
function fftangle(img; title="phase")
    heatmap(angle.(fftshift(fft(img))), color=:grays, xticks=false, yticks=false, title=title)
end

"""
    fftimggui()

Build interactive GUI for playing with Fourier transform and image filters.
"""
function fftimggui()
    images = OrderedDict(
        "mandrill" => Gray.(testimage("mandrill")),
        "cameraman" => Gray.(testimage("cameraman")),
        "house" => Gray.(testimage("house"))
    )
    filters = OrderedDict(
        "Dirac" => centered([0.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0]),
        "Gaussian1" => Kernel.gaussian(1.0),
        "Gaussian10" => Kernel.gaussian(10.0),
        "DoG" => Kernel.DoG(2.0),
        "sobel_x" => Kernel.sobel()[2],
        "sobel_y" => Kernel.sobel()[1]
    )

    ui = @manipulate for image = images, filter = filters
        img = channelview(image)[end:-1:1, :]
        result = imfilter(img, filter)
        vbox(
            hbox(
                fftimshow(img),
                fftpower(img),
                fftangle(img)
            ),
            hbox(
                fftimshow(parent(filter), title="filter"),
                fftpower(parent(filter)),
                fftangle(parent(filter))
            ),
            hbox(
                fftimshow(result, title="result"),
                fftpower(result),
                fftangle(result)
            ),
        )
    end
    return ui
end

end # module Utils
