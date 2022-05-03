#include "AudioThumbnailBars.h"



    struct AudioThumbnailBars::MinMaxValue
    {
        MinMaxValue() noexcept
        {
            values[0] = 0;
            values[1] = 0;
        }

        inline void set (const juce::int8 newMin, const juce::int8 newMax) noexcept
        {
            values[0] = newMin;
            values[1] = newMax;
        }

        inline juce::int8 getMinValue() const noexcept        { return values[0]; }
        inline juce::int8 getMaxValue() const noexcept        { return values[1]; }

        inline void setFloat (juce::Range<float> newRange) noexcept
        {
            // Workaround for an ndk armeabi compiler bug which crashes on signed saturation
#if JUCE_ANDROID
            Range<float> limitedRange (jlimit (-1.0f, 1.0f, newRange.getStart()),
                                   jlimit (-1.0f, 1.0f, newRange.getEnd()));
        values[0] = (int8) (limitedRange.getStart() * 127.0f);
        values[1] = (int8) (limitedRange.getEnd()   * 127.0f);
#else
            values[0] = (juce::int8) juce::jlimit (-128, 127, juce::roundToInt (newRange.getStart() * 127.0f));
            values[1] = (juce::int8) juce::jlimit (-128, 127, juce::roundToInt (newRange.getEnd()   * 127.0f));
#endif

            if (values[0] == values[1])
            {
                if (values[1] == 127)
                    values[0]--;
                else
                    values[1]++;
            }
        }

        inline bool isNonZero() const noexcept
        {
            return values[1] > values[0];
        }

        inline int getPeak() const noexcept
        {
            return juce::jmax (std::abs ((int) values[0]),
                    std::abs ((int) values[1]));
        }

        inline void read (juce::InputStream& input)      { input.read (values, 2); }
        inline void write (juce::OutputStream& output)   { output.write (values, 2); }

    private:
        juce::int8 values[2];
    };

//==============================================================================
    class AudioThumbnailBars::LevelDataSource   : public juce::TimeSliceClient
    {
    public:
        LevelDataSource (AudioThumbnailBars& thumb, juce::AudioFormatReader* newReader, juce::int64 hash)
                : hashCode (hash), owner (thumb), reader (newReader)
        {
        }

        LevelDataSource (AudioThumbnailBars& thumb, juce::InputSource* src)
                : hashCode (src->hashCode()), owner (thumb), source (src)
        {
        }

        ~LevelDataSource() override
        {
            owner.cache.getTimeSliceThread().removeTimeSliceClient (this);
        }

        enum { timeBeforeDeletingReader = 3000 };

        void initialise (juce::int64 samplesFinished)
        {
            const juce::ScopedLock sl (readerLock);

            numSamplesFinished = samplesFinished;

            createReader();

            if (reader != nullptr)
            {
                lengthInSamples = reader->lengthInSamples;
                numChannels = reader->numChannels;
                sampleRate = reader->sampleRate;

                if (lengthInSamples <= 0 || isFullyLoaded())
                    reader.reset();
                else
                    owner.cache.getTimeSliceThread().addTimeSliceClient (this);
            }
        }

        void getLevels (juce::int64 startSample, int numSamples, juce::Array<juce::Range<float>>& levels)
        {
            const juce::ScopedLock sl (readerLock);

            if (reader == nullptr)
            {
                createReader();

                if (reader != nullptr)
                {
                    lastReaderUseTime = juce::Time::getMillisecondCounter();
                    owner.cache.getTimeSliceThread().addTimeSliceClient (this);
                }
            }

            if (reader != nullptr)
            {
                if (levels.size() < (int) reader->numChannels)
                    levels.insertMultiple (0, {}, (int) reader->numChannels - levels.size());

                reader->readMaxLevels (startSample, numSamples, levels.getRawDataPointer(), (int) reader->numChannels);

                lastReaderUseTime = juce::Time::getMillisecondCounter();
            }
        }

        void releaseResources()
        {
            const juce::ScopedLock sl (readerLock);
            reader.reset();
        }

        int useTimeSlice() override
        {
            if (isFullyLoaded())
            {
                if (reader != nullptr && source != nullptr)
                {
                    if (juce::Time::getMillisecondCounter() > lastReaderUseTime + timeBeforeDeletingReader)
                        releaseResources();
                    else
                        return 200;
                }

                return -1;
            }

            bool justFinished = false;

            {
                const juce::ScopedLock sl (readerLock);
                createReader();

                if (reader != nullptr)
                {
                    if (! readNextBlock())
                        return 0;

                    justFinished = true;
                }
            }

            if (justFinished)
                owner.cache.storeThumb (owner, hashCode);

            return 200;
        }

        bool isFullyLoaded() const noexcept
        {
            return numSamplesFinished >= lengthInSamples;
        }

        inline int sampleToThumbSample (const juce::int64 originalSample) const noexcept
        {
            return (int) (originalSample / owner.samplesPerThumbSample);
        }

        juce::int64 lengthInSamples = 0, numSamplesFinished = 0;
        double sampleRate = 0;
        unsigned int numChannels = 0;
        juce::int64 hashCode = 0;

    private:
        AudioThumbnailBars& owner;
        std::unique_ptr<juce::InputSource> source;
        std::unique_ptr<juce::AudioFormatReader> reader;
        juce::CriticalSection readerLock;
        std::atomic<juce::uint32> lastReaderUseTime { 0 };

        void createReader()
        {
            if (reader == nullptr && source != nullptr)
                if (auto* audioFileStream = source->createInputStream())
                    reader.reset (owner.formatManagerToUse.createReaderFor (std::unique_ptr<juce::InputStream> (audioFileStream)));
        }

        bool readNextBlock()
        {
            jassert (reader != nullptr);

            if (! isFullyLoaded())
            {
                auto numToDo = (int) juce::jmin (256 * (juce::int64) owner.samplesPerThumbSample, lengthInSamples - numSamplesFinished);

                if (numToDo > 0)
                {
                    auto startSample = numSamplesFinished;

                    auto firstThumbIndex = sampleToThumbSample (startSample);
                    auto lastThumbIndex  = sampleToThumbSample (startSample + numToDo);
                    auto numThumbSamps = lastThumbIndex - firstThumbIndex;

                    juce::HeapBlock<MinMaxValue> levelData ((unsigned int) numThumbSamps * numChannels);
                    juce::HeapBlock<MinMaxValue*> levels (numChannels);

                    for (int i = 0; i < (int) numChannels; ++i)
                        levels[i] = levelData + i * numThumbSamps;

                    juce::HeapBlock<juce::Range<float>> levelsRead (numChannels);

                    for (int i = 0; i < numThumbSamps; ++i)
                    {
                        reader->readMaxLevels ((firstThumbIndex + i) * owner.samplesPerThumbSample,
                                owner.samplesPerThumbSample, levelsRead, (int) numChannels);

                        for (int j = 0; j < (int) numChannels; ++j)
                            levels[j][i].setFloat (levelsRead[j]);
                    }

                    {
                        const juce::ScopedUnlock su (readerLock);
                        owner.setLevels (levels, firstThumbIndex, (int) numChannels, numThumbSamps);
                    }

                    numSamplesFinished += numToDo;
                    lastReaderUseTime = juce::Time::getMillisecondCounter();
                }
            }

            return isFullyLoaded();
        }
    };

//==============================================================================
    class AudioThumbnailBars::ThumbData
    {
    public:
        ThumbData (int numThumbSamples)
        {
            ensureSize (numThumbSamples);
        }

        inline MinMaxValue* getData (int thumbSampleIndex) noexcept
        {
            jassert (thumbSampleIndex < data.size());
            return data.getRawDataPointer() + thumbSampleIndex;
        }

        int getSize() const noexcept
        {
            return data.size();
        }

        void getMinMax (int startSample, int endSample, MinMaxValue& result) const noexcept
        {
            if (startSample >= 0)
            {
                endSample = juce::jmin (endSample, data.size() - 1);

                juce::int8 mx = -128;
                juce::int8 mn = 127;

                while (startSample <= endSample)
                {
                    auto& v = data.getReference (startSample);

                    if (v.getMinValue() < mn)  mn = v.getMinValue();
                    if (v.getMaxValue() > mx)  mx = v.getMaxValue();

                    ++startSample;
                }

                if (mn <= mx)
                {
                    result.set (mn, mx);
                    return;
                }
            }

            result.set (1, 0);
        }

        void write (const MinMaxValue* values, int startIndex, int numValues)
        {
            resetPeak();

            if (startIndex + numValues > data.size())
                ensureSize (startIndex + numValues);

            auto* dest = getData (startIndex);

            for (int i = 0; i < numValues; ++i)
                dest[i] = values[i];
        }

        void resetPeak() noexcept
        {
            peakLevel = -1;
        }

        int getPeak() noexcept
        {
            if (peakLevel < 0)
            {
                for (auto& s : data)
                {
                    auto peak = s.getPeak();

                    if (peak > peakLevel)
                        peakLevel = peak;
                }
            }

            return peakLevel;
        }

    private:
        juce::Array<MinMaxValue> data;
        int peakLevel = -1;

        void ensureSize (int thumbSamples)
        {
            auto extraNeeded = thumbSamples - data.size();

            if (extraNeeded > 0)
                data.insertMultiple (-1, MinMaxValue(), extraNeeded);
        }
    };

//==============================================================================
    class AudioThumbnailBars::CachedWindow
    {
    public:
        CachedWindow() {}

        void invalidate()
        {
            cacheNeedsRefilling = true;
        }

        void drawChannel (juce::Graphics& g, const juce::Rectangle<int>& area,
                const double startTime, const double endTime,
                const int channelNum, const float verticalZoomFactor,
                const double rate, const int numChans, const int sampsPerThumbSample,
                LevelDataSource* levelData, const juce::OwnedArray<ThumbData>& chans)
        {
            if (refillCache (area.getWidth(), startTime, endTime, rate,
                    numChans, sampsPerThumbSample, levelData, chans)
                    && juce::isPositiveAndBelow (channelNum, numChannelsCached))
            {
                auto clip = g.getClipBounds().getIntersection (area.withWidth (juce::jmin (numSamplesCached, area.getWidth())));

                if (! clip.isEmpty())
                {
                    auto topY = (float) area.getY();
                    auto bottomY = (float) area.getBottom();
                    auto midY = (topY + bottomY) * 0.5f;
                    auto vscale = verticalZoomFactor * (bottomY - topY) / 256.0f;

                    auto* cacheData = getData (channelNum, clip.getX() - area.getX());

                    juce::RectangleList<float> waveform;
                    waveform.ensureStorageAllocated (clip.getWidth());

                    auto x = (float) clip.getX();

                    int skipStep = 5;
                    int skipCount = 0;
                    for (int w = clip.getWidth(); --w >= 0;)
                    {
                        if (cacheData->isNonZero() && skipCount < 1)
                        {
                            auto top    = juce::jmax (midY - cacheData->getMaxValue() * vscale - 0.3f, topY);
                            auto bottom = juce::jmin (midY - cacheData->getMinValue() * vscale + 0.3f, bottomY);

                            waveform.addWithoutMerging (juce::Rectangle<float> (x, top, 1.5f, bottom - top));
                        }

                        ++skipCount;
                        skipCount = skipCount % skipStep;
                        x += 1.0f;
                        ++cacheData;
                    }

                    g.fillRectList (waveform);
                }
            }
        }

    private:
        juce::Array<MinMaxValue> data;
        double cachedStart = 0, cachedTimePerPixel = 0;
        int numChannelsCached = 0, numSamplesCached = 0;
        bool cacheNeedsRefilling = true;

        bool refillCache (int numSamples, double startTime, double endTime,
                double rate, int numChans, int sampsPerThumbSample,
                LevelDataSource* levelData, const juce::OwnedArray<ThumbData>& chans)
        {
            auto timePerPixel = (endTime - startTime) / numSamples;

            if (numSamples <= 0 || timePerPixel <= 0.0 || rate <= 0)
            {
                invalidate();
                return false;
            }

            if (numSamples == numSamplesCached
                    && numChannelsCached == numChans
                    && startTime == cachedStart
                    && timePerPixel == cachedTimePerPixel
                    && ! cacheNeedsRefilling)
            {
                return ! cacheNeedsRefilling;
            }

            numSamplesCached = numSamples;
            numChannelsCached = numChans;
            cachedStart = startTime;
            cachedTimePerPixel = timePerPixel;
            cacheNeedsRefilling = false;

            ensureSize (numSamples);

            if (timePerPixel * rate <= sampsPerThumbSample && levelData != nullptr)
            {
                auto sample = juce::roundToInt (startTime * rate);
                juce::Array<juce::Range<float>> levels;

                int i;
                for (i = 0; i < numSamples; ++i)
                {
                    auto nextSample = juce::roundToInt ((startTime + timePerPixel) * rate);

                    if (sample >= 0)
                    {
                        if (sample >= levelData->lengthInSamples)
                        {
                            for (int chan = 0; chan < numChannelsCached; ++chan)
                                *getData (chan, i) = MinMaxValue();
                        }
                        else
                        {
                            levelData->getLevels (sample, juce::jmax (1, nextSample - sample), levels);

                            auto totalChans = juce::jmin (levels.size(), numChannelsCached);

                            for (int chan = 0; chan < totalChans; ++chan)
                                getData (chan, i)->setFloat (levels.getReference (chan));
                        }
                    }

                    startTime += timePerPixel;
                    sample = nextSample;
                }

                numSamplesCached = i;
            }
            else
            {
                jassert (chans.size() == numChannelsCached);

                for (int channelNum = 0; channelNum < numChannelsCached; ++channelNum)
                {
                    ThumbData* channelData = chans.getUnchecked (channelNum);
                    MinMaxValue* cacheData = getData (channelNum, 0);

                    auto timeToThumbSampleFactor = rate / (double) sampsPerThumbSample;

                    startTime = cachedStart;
                    auto sample = juce::roundToInt (startTime * timeToThumbSampleFactor);

                    for (int i = numSamples; --i >= 0;)
                    {
                        auto nextSample = juce::roundToInt ((startTime + timePerPixel) * timeToThumbSampleFactor);

                        channelData->getMinMax (sample, nextSample, *cacheData);

                        ++cacheData;
                        startTime += timePerPixel;
                        sample = nextSample;
                    }
                }
            }

            return true;
        }

        MinMaxValue* getData (const int channelNum, const int cacheIndex) noexcept
        {
            jassert (juce::isPositiveAndBelow (channelNum, numChannelsCached) && juce::isPositiveAndBelow (cacheIndex, data.size()));

            return data.getRawDataPointer() + channelNum * numSamplesCached
                    + cacheIndex;
        }

        void ensureSize (const int numSamples)
        {
            auto itemsRequired = numSamples * numChannelsCached;

            if (data.size() < itemsRequired)
                data.insertMultiple (-1, MinMaxValue(), itemsRequired - data.size());
        }
    };

//==============================================================================
    AudioThumbnailBars::AudioThumbnailBars (const int originalSamplesPerThumbnailSample,
            juce::AudioFormatManager& formatManager,
            juce::AudioThumbnailCache& cacheToUse)
            : formatManagerToUse (formatManager),
              cache (cacheToUse),
              window (new CachedWindow()),
              samplesPerThumbSample (originalSamplesPerThumbnailSample)
    {
    }

    AudioThumbnailBars::~AudioThumbnailBars()
    {
        clear();
    }

    void AudioThumbnailBars::clear()
    {
        source.reset();
        const juce::ScopedLock sl (lock);
        clearChannelData();
    }

    void AudioThumbnailBars::clearChannelData()
    {
        window->invalidate();
        channels.clear();
        totalSamples = numSamplesFinished = 0;
        numChannels = 0;
        sampleRate = 0;

        sendChangeMessage();
    }

    void AudioThumbnailBars::reset (int newNumChannels, double newSampleRate, juce::int64 totalSamplesInSource)
    {
        clear();

        const juce::ScopedLock sl (lock);
        numChannels = newNumChannels;
        sampleRate = newSampleRate;
        totalSamples = totalSamplesInSource;

        createChannels (1 + (int) (totalSamplesInSource / samplesPerThumbSample));
    }

    void AudioThumbnailBars::createChannels (const int length)
    {
        while (channels.size() < numChannels)
            channels.add (new ThumbData (length));
    }

//==============================================================================
    bool AudioThumbnailBars::loadFrom (juce::InputStream& rawInput)
    {
        juce::BufferedInputStream input (rawInput, 4096);

        if (input.readByte() != 'j' || input.readByte() != 'a' || input.readByte() != 't' || input.readByte() != 'm')
            return false;

        const juce::ScopedLock sl (lock);
        clearChannelData();

        samplesPerThumbSample = input.readInt();
        totalSamples = input.readInt64();             // Total number of source samples.
        numSamplesFinished = input.readInt64();       // Number of valid source samples that have been read into the thumbnail.
        juce::int32 numThumbnailSamples = input.readInt();  // Number of samples in the thumbnail data.
        numChannels = input.readInt();                // Number of audio channels.
        sampleRate = input.readInt();                 // Source sample rate.
        input.skipNextBytes (16);                     // (reserved)

        createChannels (numThumbnailSamples);

        for (int i = 0; i < numThumbnailSamples; ++i)
            for (int chan = 0; chan < numChannels; ++chan)
                channels.getUnchecked(chan)->getData(i)->read (input);

        return true;
    }

    void AudioThumbnailBars::saveTo (juce::OutputStream& output) const
    {
        const juce::ScopedLock sl (lock);

        const int numThumbnailSamples = channels.size() == 0 ? 0 : channels.getUnchecked(0)->getSize();

        output.write ("jatm", 4);
        output.writeInt (samplesPerThumbSample);
        output.writeInt64 (totalSamples);
        output.writeInt64 (numSamplesFinished);
        output.writeInt (numThumbnailSamples);
        output.writeInt (numChannels);
        output.writeInt ((int) sampleRate);
        output.writeInt64 (0);
        output.writeInt64 (0);

        for (int i = 0; i < numThumbnailSamples; ++i)
            for (int chan = 0; chan < numChannels; ++chan)
                channels.getUnchecked(chan)->getData(i)->write (output);
    }

//==============================================================================
    bool AudioThumbnailBars::setDataSource (LevelDataSource* newSource)
    {
        JUCE_ASSERT_MESSAGE_MANAGER_IS_LOCKED

                numSamplesFinished = 0;
        auto wasSuccessful = [&] { return sampleRate > 0 && totalSamples > 0; };

        if (cache.loadThumb (*this, newSource->hashCode) && isFullyLoaded())
        {
            source.reset (newSource); // (make sure this isn't done before loadThumb is called)

            source->lengthInSamples = totalSamples;
            source->sampleRate = sampleRate;
            source->numChannels = (unsigned int) numChannels;
            source->numSamplesFinished = numSamplesFinished;

            return wasSuccessful();
        }

        source.reset (newSource);

        const juce::ScopedLock sl (lock);
        source->initialise (numSamplesFinished);

        totalSamples = source->lengthInSamples;
        sampleRate = source->sampleRate;
        numChannels = (juce::int32) source->numChannels;

        createChannels (1 + (int) (totalSamples / samplesPerThumbSample));

        return wasSuccessful();
    }

    bool AudioThumbnailBars::setSource (juce::InputSource* const newSource)
    {
        clear();

        return newSource != nullptr && setDataSource (new LevelDataSource (*this, newSource));
    }

    void AudioThumbnailBars::setReader (juce::AudioFormatReader* newReader, juce::int64 hash)
    {
        clear();

        if (newReader != nullptr)
            setDataSource (new LevelDataSource (*this, newReader, hash));
    }

juce::int64 AudioThumbnailBars::getHashCode() const
    {
        return source == nullptr ? 0 : source->hashCode;
    }

    void AudioThumbnailBars::addBlock (juce::int64 startSample, const juce::AudioBuffer<float>& incoming,
            int startOffsetInBuffer, int numSamples)
    {
        jassert (startSample >= 0
                && startOffsetInBuffer >= 0
                && startOffsetInBuffer + numSamples <= incoming.getNumSamples());

        auto firstThumbIndex = (int) (startSample / samplesPerThumbSample);
        auto lastThumbIndex  = (int) ((startSample + numSamples + (samplesPerThumbSample - 1)) / samplesPerThumbSample);
        auto numToDo = lastThumbIndex - firstThumbIndex;

        if (numToDo > 0)
        {
            auto numChans = juce::jmin (channels.size(), incoming.getNumChannels());

            const juce::HeapBlock<MinMaxValue> thumbData (numToDo * numChans);
            const juce::HeapBlock<MinMaxValue*> thumbChannels (numChans);

            for (int chan = 0; chan < numChans; ++chan)
            {
                auto* sourceData = incoming.getReadPointer (chan, startOffsetInBuffer);
                auto* dest = thumbData + numToDo * chan;
                thumbChannels [chan] = dest;

                for (int i = 0; i < numToDo; ++i)
                {
                    auto start = i * samplesPerThumbSample;
                    dest[i].setFloat (juce::FloatVectorOperations::findMinAndMax (sourceData + start, juce::jmin (samplesPerThumbSample, numSamples - start)));
                }
            }

            setLevels (thumbChannels, firstThumbIndex, numChans, numToDo);
        }
    }

    void AudioThumbnailBars::setLevels (const MinMaxValue* const* values, int thumbIndex, int numChans, int numValues)
    {
        const juce::ScopedLock sl (lock);

        for (int i = juce::jmin (numChans, channels.size()); --i >= 0;)
            channels.getUnchecked(i)->write (values[i], thumbIndex, numValues);

        auto start = thumbIndex * (juce::int64) samplesPerThumbSample;
        auto end   = (thumbIndex + numValues) * (juce::int64) samplesPerThumbSample;

        if (numSamplesFinished >= start && end > numSamplesFinished)
            numSamplesFinished = end;

        totalSamples = juce::jmax (numSamplesFinished, totalSamples.load());
        window->invalidate();
        sendChangeMessage();
    }

//==============================================================================
    int AudioThumbnailBars::getNumChannels() const noexcept
    {
        return numChannels;
    }

    double AudioThumbnailBars::getTotalLength() const noexcept
    {
        return sampleRate > 0 ? ((double) totalSamples / sampleRate) : 0.0;
    }

    bool AudioThumbnailBars::isFullyLoaded() const noexcept
    {
        return numSamplesFinished >= totalSamples - samplesPerThumbSample;
    }

    double AudioThumbnailBars::getProportionComplete() const noexcept
    {
        return juce::jlimit (0.0, 1.0, (double) numSamplesFinished / (double) juce::jmax ((juce::int64) 1, totalSamples.load()));
    }

juce::int64 AudioThumbnailBars::getNumSamplesFinished() const noexcept
    {
        return numSamplesFinished;
    }

    float AudioThumbnailBars::getApproximatePeak() const
    {
        const juce::ScopedLock sl (lock);
        int peak = 0;

        for (auto* c : channels)
            peak = juce::jmax (peak, c->getPeak());

        return (float) juce::jlimit (0, 127, peak) / 127.0f;
    }

    void AudioThumbnailBars::getApproximateMinMax (double startTime, double endTime, int channelIndex,
            float& minValue, float& maxValue) const noexcept
    {
        const juce::ScopedLock sl (lock);
        MinMaxValue result;
        auto* data = channels [channelIndex];

        if (data != nullptr && sampleRate > 0)
        {
            auto firstThumbIndex = (int) ((startTime * sampleRate) / samplesPerThumbSample);
            auto lastThumbIndex  = (int) (((endTime * sampleRate) + samplesPerThumbSample - 1) / samplesPerThumbSample);

            data->getMinMax (juce::jmax (0, firstThumbIndex), lastThumbIndex, result);
        }

        minValue = result.getMinValue() / 128.0f;
        maxValue = result.getMaxValue() / 128.0f;
    }

    void AudioThumbnailBars::drawChannel (juce::Graphics& g, const juce::Rectangle<int>& area, double startTime,
            double endTime, int channelNum, float verticalZoomFactor)
    {
        const juce::ScopedLock sl (lock);

        window->drawChannel (g, area, startTime, endTime, channelNum, verticalZoomFactor,
                sampleRate, numChannels, samplesPerThumbSample, source.get(), channels);
    }

    void AudioThumbnailBars::drawChannels (juce::Graphics& g, const juce::Rectangle<int>& area, double startTimeSeconds,
            double endTimeSeconds, float verticalZoomFactor)
    {
        for (int i = 0; i < numChannels; ++i)
        {
            auto y1 = juce::roundToInt ((i * area.getHeight()) / numChannels);
            auto y2 = juce::roundToInt (((i + 1) * area.getHeight()) / numChannels);

            drawChannel (g, { area.getX(), area.getY() + y1, area.getWidth(), y2 - y1 },
                    startTimeSeconds, endTimeSeconds, i, verticalZoomFactor);
        }
    }
